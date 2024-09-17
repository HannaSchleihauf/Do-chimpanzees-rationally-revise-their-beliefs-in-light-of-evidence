# Study chimpanzee: Do chimpanzees rationally revise their beliefs in light of
# evidence?
# Analyses Study 4 - Second-Order Evidence

# Last modified: Sep 10, 2024 5:55

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
# load data
dat4 <-
  read.csv(
    file =
      "./study4/data/data_study_4.csv",
    head = TRUE, stringsAsFactors = TRUE
  )

############################################################################
# DATA WRANGELING
############################################################################
dat4 <-
  dat4 %>%
  filter(!grepl("filler", condition)) # deleted fillers
dat4 <-
  dat4 %>%
  arrange(chimpanzee, trial) %>%
  group_by(condition) %>%
  mutate(trial.per.condition = row_number())

############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
xx.fe.re <- fe.re.tab(
  fe.model =
    "belief.revision ~ condition*trial.per.condition",
  re = "(1|chimpanzee)",
  data = dat4
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

# Center dummy variables
t.data$condition.independent.code <-
  t.data$condition.independent -
  mean(t.data$condition.independent)
t.data$z.trial <-
  scale(t.data$trial.per.condition)

############################################################################
# FITTING THE MODEL
############################################################################
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000000))

full_exp4 <- glmer(
  belief.revision ~
    condition * z.trial +
    (1 + condition.independent.code * z.trial || chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)

main_exp4 <- glmer(
  belief.revision ~
    condition + z.trial +
    (1 + condition.independent.code * z.trial | chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)

null_exp4 <- glmer(
  belief.revision ~
    1 +
    (1 + condition.independent.code * z.trial | chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)

summary(full_exp4)$varcor
ranef.diagn.plot(full_exp4)
round(summary(full_exp4)$coefficients, 3)

############################################################################
# CHECKING ASSUMPTIONS
############################################################################
overdisp.test(full_exp4)
library(car)
vif(main_exp4)
# Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full_exp4, contr = contr, use = c("chimpanzee"))
m.stab.b$detailed$warnings
as.data.frame(round(m.stab.b$summary[, -1], 3))
m.stab.plot(round(m.stab.b$summary[, -1], 3))

############################################################################
# MODEL COMPARISONS
############################################################################
round(anova(full_exp4, null_exp4, test = "Chisq"), 3)
round(drop1(full_exp4, test = "Chisq"), 3)
round(drop1(main_exp4, test = "Chisq"), 3)

## First peek at effects
plot(effect("condition*z.trial", full_exp4), type = "response")
plot(effect("condition", main_exp4), type = "response")
plot(effect("z.trial", main_exp4), type = "response")

# Coefficients of the full_exp4 model
round(summary(full_exp4)$coefficients, 3)

############################################################################
# BOOTSTRAPS
############################################################################
## Bootstraps of full_exp4 model
# The bootstrap has already been run and is saved in the image
boot.full_exp4 <- boot.glmm.pred(
  model.res = full_exp4, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("condition", "z.trial")
)
round(boot.full_exp4$ci.estimates, 3)
m.stab.plot(round(boot.full_exp4$ci.estimates, 3))
boot.full_exp4$ci.predicted

boot.main_exp4 <- boot.glmm.pred(
  model.res = main_exp4, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("condition")
)
round(boot.main_exp4$ci.estimates, 3)
as.data.frame(round(boot.main_exp4$ci.estimates, 3))
m.stab.plot(round(boot.main_exp4$ci.estimates, 3))
boot.main_exp4$ci.predicted

save.image("./R_images/study_4_analysis.RData")

############################################################################
# PLOTTING
############################################################################
load("./R_images/study_4_analysis.RData")

# Condition ONLY ###########################################################
dat4.agg <- dat4 %>%
  mutate(belief.revision.numeric = as.numeric(belief.revision) - 1) %>%
  group_by(chimpanzee, condition) %>%
  summarise(mean.resp = mean(belief.revision.numeric, na.rm = T)) %>%
  ungroup()
dat4.agg$condition <- droplevels(dat4.agg$condition)

# Data manipulation outside ggplot
dat4.agg$condition2 <-
  jitter(as.numeric(as.factor(dat4.agg$condition)), amount = 0.12)
dat4.agg$mean.resp2 <-
  jitter(dat4.agg$mean.resp, amount = 0.04)

ci_predicted_dependent <- boot.main_exp4$ci.predicted %>%
  filter(condition == "dependent")
ci_predicted_independent <- boot.main_exp4$ci.predicted %>%
  filter(condition == "independent")

# ggplot
exp4_plot_condition <-
  ggplot() +
  geom_point(
    data = dat4.agg,
    aes(x = condition2, y = mean.resp2, color = condition),
    size = 1.75, alpha = .4
  ) +
  scale_color_manual(values =
                       c("dependent" = "#4291F8", "independent" = "#F19134")) +
  geom_line(
    data = dat4.agg,
    aes(x = condition2, y = mean.resp2, group = chimpanzee),
    color = "gray", lty = 1, alpha = .7
  ) +
  geom_errorbar(
    data = ci_predicted_dependent,
    aes(
      x = as.numeric(condition) - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = condition, color = condition
    ),
    width = 0.15, #size = 1
  ) +
  geom_errorbar(
    data = ci_predicted_independent,
    aes(
      x = as.numeric(condition) - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = condition, color = condition
    ),
    width = 0.15, #size = 1
  ) +
  geom_point(
    data = ci_predicted_dependent,
    aes(x = as.numeric(condition) - 0.25, y = fitted),
    color = "#4291F8", size = 2
  ) +
  geom_point(
    data = ci_predicted_independent,
    aes(x = as.numeric(condition) - 0.25, y = fitted),
    color = "#F19134", size = 2
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "gray"
  ) + # Adding the dashed line here
  scale_x_discrete(
    limits = c("dependent", "independent"),
    name = "evidence condition",
    labels = c("dependent", "independent")
  ) +
  scale_y_continuous(
    name = "proportion of belief revision",
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    limits = c(-0.05, 1),
    labels = scales::percent
  ) +
  geom_violin(
    data = dat4.agg,
    aes(x = condition, y = mean.resp, fill = condition),
    position = position_nudge(x = 0),
    alpha = .2
  ) +
  scale_fill_manual(values = c(
    "dependent" = "#4291F8",
    "independent" = "#F19134"
  )) +
  labs(title = "Experiment 4") +
  theme_classic(base_size=12)+
  theme(legend.position = "none",
        strip.text.x = element_text(hjust = 0),
        axis.title.y = element_text(margin=margin(t=0,r=19,b=0,l=0)),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.y.left = element_text(vjust = -2.4),
        axis.title.x = element_text(
          margin =
            margin(t = 10, r = 0, b = 0, l = 0)
        ))

exp4_plot_condition

save.image("./R_images/study_4_analysis.RData")
