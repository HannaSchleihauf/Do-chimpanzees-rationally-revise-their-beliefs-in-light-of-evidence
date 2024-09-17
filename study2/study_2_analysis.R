# Study Name: Do chimpanzees rationally revise their beliefs in light of evidence?
# Stduy 1

# Last modified: Jan 8, 2024 21:17

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
# LOAD DATA & DATA WRANGELING
############################################################################
xdata <-
  read.csv(
    file =
      "./study2/data/Results_Belief_Revision.csv",
    head = TRUE, stringsAsFactors = TRUE
  )
str(xdata)
xdata$Belief_Revision[xdata$Belief_Revision == "/"] <- NA
xdata$Belief_Revision <- droplevels(xdata$Belief_Revision)

xdata <-
  xdata %>%
  filter(!grepl("filler", Condition)) # deleted fillers
xdata <-
  xdata %>%
  arrange(Chimpanzee, Trial) %>%
  group_by(Chimpanzee, Condition) %>%
  mutate(trial.per.Condition = row_number()) %>%
  ungroup()

tapply(as.numeric(xdata$Belief_Revision)-1,
       list(xdata$Condition, as.factor(xdata$trial.per.Condition)), mean)


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
contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))

contr = glmerControl(optimizer = "nloptwrap",
optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD",  maxit = 1e9))

full_exp2 <- glmer(Belief_Revision ~
                Condition * z.Trial +
                (1 + Condition.weak_first.code * z.Trial | Chimpanzee),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

main_exp2 <- glmer(Belief_Revision ~
                Condition + z.Trial +
                (1 + Condition.weak_first.code * z.Trial | Chimpanzee),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

null_exp2 <- glmer(Belief_Revision ~
                1 +
                (1 + Condition.weak_first.code * z.Trial | Chimpanzee),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

summary(full_exp2)$varcor
ranef.diagn.plot(full_exp2)
round(summary(full_exp2)$coefficients, 3)

############################################################################
# CHECKING ASSUMPTIONS
############################################################################
overdisp.test(full_exp2)
library(car)
vif(main_exp2)
# Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full_exp2, contr = contr, use = c("Chimpanzee"))
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
round(anova(full_exp2, null_exp2, test = "Chisq"), 3)
round(drop1(full_exp2, test = "Chisq"), 3)
round(drop1(main_exp2, test = "Chisq"), 3)

## First peek at effects
plot(effect("Condition", main_exp2), type = "response")
plot(effect("z.Trial", main_exp2), type = "response")

# Coefficients of the full_exp2 model
round(summary(full_exp2)$coefficients, 3)

############################################################################
# BOOTSTRAPS
############################################################################
## Bootstraps of full_exp2 model
# The bootstrap has already been run and is saved in the image
boot.full_exp2 <- boot.glmm.pred(
  model.res = full_exp2, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("Condition", "z.Trial")
)

round(boot.full_exp2$ci.estimates, 3)
m.stab.plot(round(boot.full_exp2$ci.estimates, 3))
boot.full_exp2$ci.predicted

boot.main_exp2 <- boot.glmm.pred(
  model.res = main_exp2, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("Condition")
)

round(boot.main_exp2$ci.estimates, 3)
m.stab.plot(round(boot.main_exp2$ci.estimates, 3))
boot.main_exp2$ci.predicted

save.image("./R_images/study_2_analysis.RData")

############################################################################
# PLOTTING
############################################################################
load("./R_images/study_2_analysis.RData")

library(gghalves)
library(ggthemes)
library(cowplot)

xdata.agg <- xdata %>%
  mutate(belief_revision.numeric = as.numeric(Belief_Revision) - 1) %>%
  group_by(Chimpanzee, Condition) %>%
  summarise(mean.resp = mean(belief_revision.numeric, na.rm = T)) %>%
  ungroup()
xdata.agg$Condition <- droplevels(xdata.agg$Condition)

# Data manipulation outside ggplot
xdata.agg$Condition2 <-
  jitter(as.numeric(as.factor(xdata.agg$Condition)), amount = 0.12)
xdata.agg$mean.resp2 <-
  jitter(xdata.agg$mean.resp, amount = 0.04)

ci_predicted_strong <- boot.main_exp2$ci.predicted %>%
  filter(Condition == "strong_first")
ci_predicted_weak <- boot.main_exp2$ci.predicted %>%
  filter(Condition == "weak_first")

# ggplot
exp2_plot_Condition <-
  ggplot() +
  geom_point(
    data = xdata.agg,
    aes(x = Condition2, y = mean.resp2, color = Condition),
    size = 2.5, alpha = .4
  ) +
  scale_color_manual(values =
    c("strong_first" = "dodgerblue", "weak_first" = "darkorange")) +
  geom_line(
    data = xdata.agg,
    aes(x = Condition2, y = mean.resp2, group = Chimpanzee),
    color = "gray", lty = 1, alpha = .7
  ) +
  geom_errorbar(
    data = ci_predicted_strong,
    aes(
      x = as.numeric(Condition) - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = Condition, color = Condition
    ),
    width = 0.1, size = 1
  ) +
  geom_errorbar(
    data = ci_predicted_weak,
    aes(
      x = as.numeric(Condition) + 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = Condition, color = Condition
    ),
    width = 0.1, size = 1
  ) +
  geom_point(
    data = ci_predicted_strong,
    aes(x = as.numeric(Condition) - 0.25, y = fitted),
    color = "dodgerblue", size = 2.5
  ) +
  geom_point(
    data = ci_predicted_weak,
    aes(x = as.numeric(Condition) + 0.25, y = fitted),
    color = "darkorange", size = 2.5
  ) +
  scale_x_discrete(
    limits = c("strong_first", "weak_first"),
    name = "Condition",
    labels = c("strong evidence first", "weak evidence first")
  ) +
  scale_y_continuous(
    name = "proportion of belief revision",
    labels = scales::percent
  ) +
  geom_violin(
    data = xdata.agg,
    aes(x = Condition, y = mean.resp, fill = Condition),
    position = position_nudge(x = 0),
    alpha = .2
  ) +
  scale_fill_manual(values = c(
    "strong_first" = "dodgerblue",
    "weak_first" = "darkorange"
  )) +
  labs(title = "Experiment 2") +
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
exp2_plot_Condition

save.image("./R_images/study_2_analysis.RData")
