# Study Chimpanzee: Do Chimpanzees rationally revise their beliefs in light of
# evidence?
# Study 1

# Last modified: Feb 27, 2024 5:55

############################################################################
# PACKAGES & FUNCTIONS
############################################################################
library("easypackages")
libraries(
  "lme4", "tidyverse", "tidyselect",
  "optimx", "emmeans", "car", "effects"
)

source("./functions/diagnostic_fcns.r")
source("./functions/glmm_stability.r")
source("./functions/drop1_para.r")
source("./functions/boot_glmm.r")

############################################################################
# R SETTINGS
############################################################################
options(scipen = 999)

############################################################################
# DATA
############################################################################
xdata <-
  read.csv(
    file =
      "./study1/data/Results_Belief_Revision.csv",
    head = TRUE, stringsAsFactors = TRUE
  )
str(xdata)
xdata$Belief_Revision[xdata$Belief_Revision == "/"] <- NA
xdata$Belief_Revision <- droplevels(xdata$Belief_Revision)

############################################################################
# DATA WRANGELING
############################################################################
xdata <-
  xdata %>%
  filter(!grepl("filler", Condition)) # deleted fillers
xdata <-
  xdata %>%
  arrange(Chimpanzee, Trial) %>%
  group_by(Condition) %>%
  mutate(trial.per.Condition = row_number())

############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
xx.fe.re <- fe.re.tab(
  fe.model =
    "Belief_Revision ~ Condition*trial.per.Condition",
  re = "(1|Chimpanzee)",
  data = xdata
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

# Center dummy variables
t.data$Condition.weak_first.code <-
  t.data$Condition.weak_first - mean(t.data$Condition.weak_first)
t.data$z.Trial <-
  scale(t.data$trial.per.Condition)

############################################################################
# FITTING THE MODEL AS PREREGISTERED
############################################################################
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000000))

full_exp1 <- glmer(Belief_Revision ~
  Condition * z.Trial +
  (1 + Condition.weak_first.code * z.Trial | Chimpanzee),
data = t.data, control = contr,
family = binomial(link = "logit")
)

main_exp1 <- glmer(Belief_Revision ~
  Condition + z.Trial +
  (1 + Condition.weak_first.code * z.Trial | Chimpanzee),
data = t.data, control = contr,
family = binomial(link = "logit")
)

null_exp1 <- glmer(Belief_Revision ~
  1 +
  (1 + Condition.weak_first.code * z.Trial | Chimpanzee),
data = t.data, control = contr,
family = binomial(link = "logit")
)

summary(full_exp1)$varcor
ranef.diagn.plot(full_exp1)
round(summary(full_exp1)$coefficients, 3)

############################################################################
# CHECKING ASSUMPTIONS
############################################################################
overdisp.test(full_exp1)
library(car)
vif(main_exp1)
# Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full_exp1, contr = contr, use = c("Chimpanzee"))
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
round(anova(full_exp1, null_exp1, test = "Chisq"), 3)
round(drop1(full_exp1, test = "Chisq"), 3)
round(drop1(main_exp1, test = "Chisq"), 3)

## First peek at effects
plot(effect("Condition", main_exp1), type = "response")
plot(effect("z.Trial", main_exp1), type = "response")

# Coefficients of the full_exp1 model
round(summary(full_exp1)$coefficients, 3)

############################################################################
# BOOTSTRAPS
############################################################################
## Bootstraps of full_exp1 model
# The bootstrap has already been run and is saved in the image
boot.full_exp1 <- boot.glmm.pred(
  model.res = full_exp1, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("Condition", "z.Trial")
)

round(boot.full_exp1$ci.estimates, 3)
as.data.frame(round(boot.main_exp1$ci.estimates, 3))
m.stab.plot(round(boot.main_exp1$ci.estimates, 3))
boot.main_exp1$ci.predicted

boot.main_exp1 <- boot.glmm.pred(
  model.res = main_exp1, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("Condition")
)

round(boot.main_exp1$ci.estimates, 3)
as.data.frame(round(boot.main_exp1$ci.estimates, 3))
m.stab.plot(round(boot.main_exp1$ci.estimates, 3))
boot.main_exp1$ci.predicted

save.image("./R_images/study_1_analysis.RData")

############################################################################
# PLOTTING
############################################################################
load("./R_images/study_1_analysis.RData")

xdata.agg <- xdata2 %>%
  mutate(first.evidence.chosen.numeric =
    as.numeric(as.factor(first.evidence.chosen)) - 1) %>%
  group_by(Chimpanzee, first.evidence) %>%
  summarise(mean.resp = mean(first.evidence.chosen.numeric, na.rm = T)) %>%
  ungroup()
xdata.agg <- droplevels(xdata.agg)

xdata2$first.evidence <-
  factor(xdata2$first.evidence,
         levels = c(
           "auditory_strong",
           "trace_weak"
         )
  )
xdata.agg$first.evidence <-
  factor(xdata.agg$first.evidence,
         levels = c(
           "auditory_strong",
           "trace_weak"
         )
  )
levels(xdata.agg$first.evidence)

# Data manipulation outside ggplot
xdata.agg$first.evidence2 <-
  jitter(as.numeric(as.factor(xdata.agg$first.evidence)), amount = 0.12)
xdata.agg$mean.resp2 <-
  jitter(xdata.agg$mean.resp, amount = 0.04)

ci_predicted_study2_trace <- boot.main_exp.2$ci.predicted %>%
  filter(first.evidence == "trace_weak")
ci_predicted_study2_auditory <- boot.main_exp.2$ci.predicted %>%
  filter(first.evidence == "auditory_strong")

# ggplot
exp2_plot_first_choice <-
  ggplot() +
  geom_point(
    data = xdata.agg,
    aes(x = first.evidence2, y = mean.resp2, color = first.evidence),
    size = 2.5, alpha = .4
  ) +
  scale_color_manual(values = c(
    "auditory_strong" = "dodgerblue",
    "trace_weak" = "darkorange"
  )) +
  geom_line(
    data = xdata.agg,
    aes(x = first.evidence2, y = mean.resp2, group = Chimpanzee),
    color = "gray", lty = 1, alpha = .7
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray30") +
  geom_errorbar(
    data = ci_predicted_study2_auditory,
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
      x = 2 - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = first.evidence, color = first.evidence
    ),
    width = 0.1, size = 1
  ) +
  geom_point(
    data = ci_predicted_study1_visual,
    aes(x = 1 - 0.25, y = fitted),
    color = "dodgerblue", size = 2.5
  ) +
  geom_point(
    data = ci_predicted_study1_auditory,
    aes(x = 2 - 0.25, y = fitted),
    color = "darkorange", size = 2.5
  ) +
  scale_x_discrete(
    limits = c(
      "auditory_strong",
      "trace_weak"
    ),
    name = "type of evidence",
    labels = c(
      "auditory (strong)",
      "trace (weak)"
    )
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
  scale_fill_manual(values = c(
    "auditory_strong" = "dodgerblue",
    "trace_weak" = "darkorange"
  )) +
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
exp2_plot_first_choice

save.image("./R_images/study_1_analysis.RData")
