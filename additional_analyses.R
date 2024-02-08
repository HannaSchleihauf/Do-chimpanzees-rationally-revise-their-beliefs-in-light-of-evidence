library("xlsx")
library("tidyverse")
library("emmeans")
library("effects")
library("lme4")

source("./functions/diagnostic_fcns.r")
source("./functions/glmm_stability.r")
source("./functions/drop1_para.r")
source("./functions/boot_glmm.r")

############################################################################
# DATA
############################################################################
xdata1 <-
  read.csv(
    file =
      "./study1/data/Results_Belief_Revision.csv",
    head = TRUE, stringsAsFactors = TRUE
  )
xdata1$Belief_Revision[xdata1$Belief_Revision == "/"] <- NA
xdata1$Belief_Revision <- droplevels(xdata1$Belief_Revision)

table(xdata1$Condition, xdata1$First_Choice)
table(xdata2$Condition, xdata2$First_Choice)

xdata1 <-
  xdata1 %>%
  filter(!grepl("filler", Condition)) %>%
  arrange(Chimpanzee, Trial) %>%
  group_by(Condition) %>%
  mutate(trial.per.Condition = row_number()) %>%
  mutate(first.evidence = case_when(
    Condition == "strong_first" ~ "Study_1_visual_evidence",
    Condition == "weak_first" ~ "Study_1_auditory_evidence"
  )) %>%
  mutate(first.evidence.chosen = case_when(
    Condition == "strong_first" & First_Choice == "strong" ~ "yes",
    Condition == "strong_first" & First_Choice == "weak" ~ "no",
    Condition == "weak_first" & First_Choice == "weak" ~ "yes",
    Condition == "weak_first" & First_Choice == "strong" ~ "no"
  ))

xdata2 <-
  read.csv(
    file =
      "./study2/data/Results_Belief_Revision.csv",
    head = TRUE, stringsAsFactors = TRUE
  )
xdata2$Belief_Revision[xdata2$Belief_Revision == "/"] <- NA
xdata2$Belief_Revision <- droplevels(xdata2$Belief_Revision)

table(xdata2$Condition, xdata2$First_Choice)

xdata2 <-
  xdata2 %>%
  filter(!grepl("filler", Condition)) %>%
  arrange(Chimpanzee, Trial) %>%
  group_by(Condition) %>%
  mutate(trial.per.Condition = row_number()) %>%
  mutate(first.evidence = case_when(
    Condition == "strong_first" ~ "Study_2_auditory_evidence",
    Condition == "weak_first" ~ "Study_2_trace_evidence"
  )) %>%
  mutate(first.evidence.chosen = case_when(
    Condition == "strong_first" & First_Choice == "strong" ~ "yes",
    Condition == "strong_first" & First_Choice == "weak" ~ "no",
    Condition == "weak_first" & First_Choice == "weak" ~ "yes",
    Condition == "weak_first" & First_Choice == "strong" ~ "no"
  ))

xdata.both <- rbind(xdata1, xdata2)

############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
xx.fe.re <- fe.re.tab(
  fe.model =
    "first.evidence.chosen ~ first.evidence + trial.per.Condition",
  re = "(1|Chimpanzee)",
  data = xdata.both
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

# Center dummy variables
t.data$first.evidence.Study_1_visual_evidence.code <-
  t.data$first.evidence.Study_1_visual_evidence -
  mean(t.data$first.evidence.Study_1_visual_evidence)
t.data$first.evidence.Study_2_auditory_evidence.code <-
  t.data$first.evidence.Study_2_auditory_evidence -
  mean(t.data$first.evidence.Study_2_auditory_evidence)
t.data$first.evidence.Study_2_trace_evidence.code <-
  t.data$first.evidence.Study_2_trace_evidence -
  mean(t.data$first.evidence.Study_2_trace_evidence)
t.data$z.Trial <-
  scale(t.data$trial.per.Condition)

############################################################################
# FITTING THE MODEL AS PREREGISTERED
############################################################################
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000000))

full_exp.both <- glmer(first.evidence.chosen ~
                     first.evidence * z.Trial +
                     (1 + (first.evidence.Study_1_visual_evidence.code +
                             first.evidence.Study_2_auditory_evidence.code +
                             first.evidence.Study_2_trace_evidence.code) *
                        z.Trial || Chimpanzee),
                   data = t.data, control = contr,
                   family = binomial(link = "logit")
)

main_exp.both <- glmer(first.evidence.chosen ~
                         first.evidence + z.Trial +
                         (1 + (first.evidence.Study_1_visual_evidence.code +
                                 first.evidence.Study_2_auditory_evidence.code +
                                 first.evidence.Study_2_trace_evidence.code) *
                            z.Trial || Chimpanzee),
                       data = t.data, control = contr,
                       family = binomial(link = "logit")
)

null_exp.both <- glmer(first.evidence.chosen ~
                         1 +
                         (1 + (first.evidence.Study_1_visual_evidence.code +
                                 first.evidence.Study_2_auditory_evidence.code +
                                 first.evidence.Study_2_trace_evidence.code) *
                            z.Trial || Chimpanzee),
                       data = t.data, control = contr,
                       family = binomial(link = "logit")
)


summary(full_exp.both)$varcor
ranef.diagn.plot(full_exp.both)
round(summary(full_exp.both)$coefficients, 3)

############################################################################
# CHECKING ASSUMPTIONS
############################################################################
overdisp.test(full_exp.both)
library(car)
vif(main_exp.both)
# Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full_exp.both, contr = contr, use = c("Chimpanzee"))
m.stab.b$detailed$warnings
as.data.frame(round(m.stab.b$summary[, -1], 3))
m.stab.plot(round(m.stab.b$summary[, -1], 3))
# Model stablility only with models that did not give a warning message
m.stab.b$detailed <- m.stab.b$detailed %>%
  filter(warnings == "none")
cbind(
  apply(m.stab.b$detailed, 2, min),
  apply(m.stab.b$detailed, 2, max)
)


############################################################################
# MODEL COMPARISONS
############################################################################
round(anova(full_exp.both, null_exp.both, test = "Chisq"), 3)
round(drop1(full_exp.both, test = "Chisq"), 3)
round(drop1(main_exp.both, test = "Chisq"), 3)

## First peek at effects
plot(effect("first.evidence", main_exp.both), type = "response")
plot(effect("z.Trial", main_exp.both), type = "response")

# Coefficients of the full_exp1 model
round(summary(full_exp.both)$coefficients, 3)

############################################################################
# BOOTSTRAPS
############################################################################

## Bootstraps of full_exp.both model
# The bootstrap has already been run and is saved in the image
boot.full_exp.both <- boot.glmm.pred(
  model.res = full_exp.both, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("first.evidence", "z.Trial")
)

round(boot.full_exp.both$ci.estimates, 3)
as.data.frame(round(boot.full_exp.both$ci.estimates, 3))
m.stab.plot(round(boot.full_exp.both$ci.estimates, 3))
boot.full_exp.both$ci.predicted

boot.main_exp.both <- boot.glmm.pred(
  model.res = main_exp.both, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("first.evidence")
)

round(boot.main_exp.both$ci.estimates, 3)
as.data.frame(round(boot.main_exp.both$ci.estimates, 3))
m.stab.plot(round(boot.main_exp.both$ci.estimates, 3))
boot.main_exp.both$ci.predicted

############################################################################
# PLOTTING
############################################################################
library(gghalves)
library(ggthemes)
library(cowplot)

xdata.agg <- xdata.both %>%
  mutate(first.evidence.chosen.numeric = as.numeric(as.factor(first.evidence.chosen)) - 1) %>%
  group_by(Chimpanzee, first.evidence) %>%
  summarise(mean.resp = mean(first.evidence.chosen.numeric, na.rm = T)) %>%
  ungroup()
xdata.agg <- droplevels(xdata.agg)

xdata.both$first.evidence <-
  factor(xdata.both$first.evidence,
         levels=c("Study_1_visual_evidence",
                  "Study_1_auditory_evidence",
                  "Study_2_auditory_evidence",
                  "Study_2_trace_evidence"))
xdata.agg$first.evidence <-
  factor(xdata.agg$first.evidence,
         levels=c( "Study_1_visual_evidence",
                   "Study_1_auditory_evidence",
                  "Study_2_auditory_evidence",
                  "Study_2_trace_evidence"))

levels(xdata.agg$first.evidence)

# Data manipulation outside ggplot
xdata.agg$first.evidence2 <-
  jitter(as.numeric(as.factor(xdata.agg$first.evidence)), amount = 0.12)

xdata.agg$mean.resp2 <-
  jitter(xdata.agg$mean.resp, amount = 0.04)



ci_predicted_study1_auditory <- boot.main_exp.both$ci.predicted %>%
  filter(first.evidence == "Study_1_auditory_evidence")
ci_predicted_study1_visual <- boot.main_exp.both$ci.predicted %>%
  filter(first.evidence == "Study_1_visual_evidence")
ci_predicted_study2_trace <- boot.main_exp.both$ci.predicted %>%
  filter(first.evidence == "Study_2_trace_evidence")
ci_predicted_study2_auditory <- boot.main_exp.both$ci.predicted %>%
  filter(first.evidence == "Study_2_auditory_evidence")

# ggplot
exp1_plot_first_choice <-
  ggplot() +
  geom_point(
    data = xdata.agg,
    aes(x = first.evidence2, y = mean.resp2, color = first.evidence),
    size = 2.5, alpha = .4
  ) +
  scale_color_manual(values = c("Study_1_visual_evidence" = "darkblue",
                                "Study_1_auditory_evidence" = "darkgreen",
                                "Study_2_auditory_evidence" = "darkviolet",
                                "Study_2_trace_evidence" = "blue")) +
  geom_line(
    data = xdata.agg,
    aes(x = first.evidence2, y = mean.resp2, group = Chimpanzee),
    color = "gray", lty = 1, alpha = .7
  ) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "gray30") +
  geom_errorbar(
    data = ci_predicted_study1_auditory,
    aes(
      x = 2 - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = first.evidence, color = first.evidence
    ),
    width = 0.1, size = 1
  ) +
  geom_errorbar(
    data = ci_predicted_study1_visual,
    aes(
      x = 1 - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = first.evidence, color = first.evidence
    ),
    width = 0.1, size = 1
  ) +
  geom_errorbar(
    data = ci_predicted_study2_trace,
    aes(
      x = 4 - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = first.evidence, color = first.evidence
    ),
    width = 0.1, size = 1
  ) +
  geom_errorbar(
    data = ci_predicted_study2_auditory,
    aes(
      x = 3 - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = first.evidence, color = first.evidence
    ),
    width = 0.1, size = 1
  ) +
  geom_point(
    data = ci_predicted_study1_visual,
    aes(x = 1 - 0.25, y = fitted),
    color = "darkblue", size = 2.5
  ) +
  geom_point(
    data = ci_predicted_study1_auditory,
    aes(x = 2 - 0.25, y = fitted),
    color = "darkgreen", size = 2.5
  ) +
  geom_point(
    data = ci_predicted_study2_auditory,
    aes(x = 3- 0.25, y = fitted),
    color = "darkviolet", size = 2.5
  ) +
  geom_point(
    data = ci_predicted_study2_trace,
    aes(x = 4 - 0.25, y = fitted),
    color = "blue", size = 2.5
  ) +
  scale_x_discrete(
    limits = c("Study_1_visual_evidence",
               "Study_1_auditory_evidence",
               "Study_2_auditory_evidence",
               "Study_2_trace_evidence"),
    name = "type of evidence",
    labels = c("Study 1\nvisual\nevidence",
               "Study 1\nauditroy\nevidence",
               "Study 2\nauditory\nevidence",
               "Study 2\ntrace\nevidence")
  ) +
  scale_y_continuous(
    name = "proportion of first choices",
    labels = scales::percent
  ) +
  geom_violin(
    data = xdata.agg,
    aes(x = first.evidence, y = mean.resp, fill = first.evidence),
    position = position_nudge(x = 0),
    alpha = .2
  ) +
  scale_fill_manual(values = c("Study_1_visual_evidence" = "darkblue",
                               "Study_1_auditory_evidence" = "darkgreen",
                               "Study_2_auditory_evidence" = "darkviolet",
                               "Study_2_trace_evidence" = "blue")) +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(
      margin =
        margin(t = 10, r = 0, b = 0, l = 0)
    ),
    legend.position = "none"
  )
exp1_plot_first_choice
