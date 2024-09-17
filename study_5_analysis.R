# Study Chimpanzee: Do Chimpanzees rationally revise their beliefs in light of
# evidence?
# Analyses Study 5 - Second-Order Evidence

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
dat5 <-
  read.csv(
    file =
      "./study5/data/data_study_5.csv",
    head = TRUE, stringsAsFactors = TRUE
  )

############################################################################
# DATA WRANGELING
############################################################################
dat5 <-
  dat5 %>%
  filter(!grepl("filler", condition)) # deleted fillers
dat5 <-
  dat5 %>%
  arrange(chimpanzee, trial) %>%
  group_by(condition) %>%
  mutate(trial.per.condition = row_number())

############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
xx.fe.re <- fe.re.tab(
  fe.model =
    "belief.revision ~ condition*trial.per.condition*domain",
  re = "(1|chimpanzee)",
  data = dat5
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

# Center dummy variables
t.data$condition.non.defeater.code <-
  t.data$condition.non.defeater -
  mean(t.data$condition.non.defeater)
t.data$z.trial <-
  scale(t.data$trial.per.condition)
t.data$domain.visual.code <-
  t.data$domain.visual -
  mean(t.data$domain.visual)

############################################################################
# FITTING THE MODEL
############################################################################
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000000))

full_exp5 <- glmer(
  belief.revision ~
    (condition * z.trial) + (condition * domain) +
    (1 + condition.non.defeater.code * z.trial +
      condition.non.defeater.code * domain.visual.code | chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)

main_exp5 <- glmer(
  belief.revision ~
    condition + z.trial + domain +
    (1 + condition.non.defeater.code * z.trial +
      condition.non.defeater.code * domain.visual.code | chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)

null_exp5 <- glmer(
  belief.revision ~
    1 +
    (1 + condition.non.defeater.code * z.trial +
      condition.non.defeater.code * domain.visual.code | chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)

summary(full_exp5)$varcor
ranef.diagn.plot(full_exp5)
round(summary(full_exp5)$coefficients, 3)

############################################################################
# CHECKING ASSUMPTIONS
############################################################################
overdisp.test(full_exp5)
library(car)
vif(main_exp5)
# Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full_exp5, contr = contr, use = c("chimpanzee"))
m.stab.b$detailed$warnings
as.data.frame(round(m.stab.b$summary[, -1], 3))
m.stab.plot(round(m.stab.b$summary[, -1], 3))

############################################################################
# MODEL COMPARISONS
############################################################################
round(anova(full_exp5, null_exp5, test = "Chisq"), 3)
round(drop1(full_exp5, test = "Chisq"), 3)
round(drop1(main_exp5, test = "Chisq"), 3)

## First peek at effects
plot(effect("condition", main_exp5), type = "response")
plot(effect("z.trial", main_exp5), type = "response")
plot(effect("domain", main_exp5), type = "response")

# Coefficients of the full_exp1 model
round(summary(full_exp5)$coefficients, 3)

############################################################################
# BOOTSTRAPS
############################################################################
## Bootstraps of full_exp5 model
# The bootstrap has already been run and is saved in the image
boot.full_exp5 <- boot.glmm.pred(
  model.res = full_exp5, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("condition", "z.trial", "domain")
)
round(boot.full_exp5$ci.estimates, 3)
m.stab.plot(round(boot.full_exp5$ci.estimates, 3))
boot.full_exp5$ci.predicted

# Bootstraps for Plots
boot.full_condition_domain <- boot.glmm.pred(
  model.res = full_exp5, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("condition", "domain")
)
round(boot.full_condition_domain$ci.estimates, 3)
m.stab.plot(round(boot.full_condition_domain$ci.estimates, 3))
boot.full_condition_domain$ci.predicted

main_exp5_condition <- glmer(
  belief.revision ~
    condition + z.trial + domain.visual.code +
    (1 + condition.non.defeater.code * z.trial +
      condition.non.defeater.code * domain.visual.code | chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)
boot.main_exp5_condition <- boot.glmm.pred(
  model.res = main_exp5_condition, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("condition")
)
round(boot.main_exp5_condition$ci.estimates, 3)
as.data.frame(round(boot.main_exp5_condition$ci.estimates, 3))
m.stab.plot(round(boot.main_exp5_condition$ci.estimates, 3))
boot.main_exp5_condition$ci.predicted

main_exp5_trial <- glmer(
  belief.revision ~
    condition.non.defeater.code + z.trial + domain.visual.code +
    (1 + condition.non.defeater.code * z.trial +
      condition.non.defeater.code * domain.visual.code | chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)
boot.main_exp5_trial <- boot.glmm.pred(
  model.res = main_exp5_trial, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("z.trial")
)
round(boot.main_exp5_trial$ci.estimates, 3)
as.data.frame(round(boot.main_exp5_trial$ci.estimates, 3))
m.stab.plot(round(boot.main_exp5_trial$ci.estimates, 3))
boot.main_exp5_trial$ci.predicted

main_exp5_domain <- glmer(
  belief.revision ~
    condition.non.defeater.code + z.trial + domain +
    (1 + condition.non.defeater.code * z.trial +
      condition.non.defeater.code * domain.visual.code | chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)
boot.main_exp5_domain <- boot.glmm.pred(
  model.res = main_exp5_domain, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("domain")
)
round(boot.main_exp5_domain$ci.estimates, 3)
as.data.frame(round(boot.main_exp5_domain$ci.estimates, 3))
m.stab.plot(round(boot.main_exp5_domain$ci.estimates, 3))
boot.main_exp5_domain$ci.predicted

save.image("./R_images/study_5_analysis.RData")

############################################################################
# PLOTTING
############################################################################
load("./R_images/study_5_analysis.RData")


# Condition AND DOMAIN #####################################################
dat5.agg_condition_domain <- dat5 %>%
  mutate(belief.revision.numeric = as.numeric(belief.revision) - 1) %>%
  group_by(chimpanzee, condition, domain) %>%
  summarise(mean.resp = mean(belief.revision.numeric, na.rm = T)) %>%
  ungroup()
dat5.agg_condition_domain$condition <- droplevels(dat5.agg$condition)

# changing factor orders
dat5.agg_condition_domain$condition <- factor(dat5.agg_condition_domain$condition,
                             levels = c("non-defeater", "defeater"))
dat5.agg_condition_domain$condition <- factor(dat5.agg_condition_domain$condition,
                         levels = c("non-defeater", "defeater"))

# Data manipulation outside ggplot
dat5.agg_condition_domain$condition2 <-
  jitter(as.numeric(as.factor(dat5.agg_condition_domain$condition)), amount = 0.12)
dat5.agg_condition_domain$mean.resp2 <-
  jitter(dat5.agg_condition_domain$mean.resp, amount = 0.04)

# changing factor orders
boot.full_condition_domain$ci.predicted$condition <-
  factor(boot.full_condition_domain$ci.predicted$condition,
         levels = c("non-defeater", "defeater"))

ci_predicted_visual <- boot.full_condition_domain$ci.predicted %>%
  filter(domain == "visual")
ci_predicted_auditory <- boot.full_condition_domain$ci.predicted %>%
  filter(domain == "auditory")

# dividing into two data sets (visual and auditory)
dat5.agg.visual <-
  subset(dat5.agg_condition_domain,
                           dat5.agg_condition_domain$domain == "visual")
dat5.agg.auditory <-
  subset(dat5.agg_condition_domain,
                             dat5.agg_condition_domain$domain == "auditory")

# ggplot
exp5_plot_visual <-
  ggplot() +
  geom_point(
    data = dat5.agg.visual,
    aes(x = condition2, y = mean.resp2, color = condition),
    size = 1.75, alpha = .4
  ) +
  scale_color_manual(values =
                       c("non-defeater" = "#4291F8", "defeater" = "#F19134")) +
  geom_line(
    data = dat5.agg.visual,
    aes(x = condition2, y = mean.resp2, group = chimpanzee),
    color = "gray", lty = 1, alpha = .7
  ) +
  geom_errorbar(
    data = ci_predicted_visual,
    aes(
      x = as.numeric(condition) - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = condition, color = condition
    ),
    width = 0.15, #size = 1
  ) +
  geom_point(
    data = ci_predicted_visual,
    aes(x = as.numeric(condition) - 0.25, y = fitted, color = condition),
    size = 2
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "gray"
  ) + # Adding the dashed line here
  scale_x_discrete(
    limits = c("non-defeater", "defeater"),
    name = "evidence condition",
    labels = c("non-defeater", "defeater")
  ) +
  scale_y_continuous(
    name = "proportion of belief revision",
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    limits = c(-0.075, 1.075),
    labels = scales::percent
  ) +
  geom_violin(
    data = dat5.agg.visual,
    aes(x = condition, y = mean.resp, fill = condition),
    position = position_nudge(x = 0),
    alpha = .2
  ) +
  scale_fill_manual(values = c(
    "non-defeater" = "#4291F8",
    "defeater" = "#F19134"
  )) +
  labs(title = "Experiment 5") +
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
exp5_plot_visual

# ggplot
exp5_plot_auditory <-
  ggplot() +
  geom_point(
    data = dat5.agg.auditory,
    aes(x = condition2, y = mean.resp2, color = condition),
    size = 1.75, alpha = .4
  ) +
  scale_color_manual(values =
                       c("non-defeater" = "#4291F8", "defeater" = "#F19134")) +
  geom_line(
    data = dat5.agg.auditory,
    aes(x = condition2, y = mean.resp2, group = chimpanzee),
    color = "gray", lty = 1, alpha = .7
  ) +
  geom_errorbar(
    data = ci_predicted_auditory,
    aes(
      x = as.numeric(condition) - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = condition, color = condition
    ),
    width = 0.15, #size = 1
  ) +
  geom_point(
    data = ci_predicted_auditory,
    aes(x = as.numeric(condition) - 0.25, y = fitted, color = condition),
    size = 2
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "gray"
  ) + # Adding the dashed line here
  scale_x_discrete(
    limits = c("non-defeater", "defeater"),
    name = "evidence condition",
    labels = c("non-defeater", "defeater")
  ) +
  scale_y_continuous(
    name = "proportion of belief revision",
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    limits = c(-0.075, 1.075),
    labels = scales::percent
  ) +
  geom_violin(
    data = dat5.agg.auditory,
    aes(x = condition, y = mean.resp, fill = condition),
    position = position_nudge(x = 0),
    alpha = .2
  ) +
  scale_fill_manual(values = c(
    "non-defeater" = "#4291F8",
    "defeater" = "#F19134"
  )) +
  labs(title = "Experiment 5") +
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
exp5_plot_auditory

# CONDITION ONLY ###########################################################
dat5.agg <- dat5 %>%
  mutate(belief.revision.numeric = as.numeric(belief.revision) - 1) %>%
  group_by(chimpanzee, condition) %>%
  summarise(mean.resp = mean(belief.revision.numeric, na.rm = T)) %>%
  ungroup()
dat5.agg$condition <- droplevels(dat5.agg$condition)
dat5.agg$condition <- factor(dat5.agg$condition,
                             levels = c("non-defeater", "defeater"))
dat5$condition <- factor(dat5$condition,
                             levels = c("non-defeater", "defeater"))

# Data manipulation outside ggplot
dat5.agg$condition2 <-
  jitter(as.numeric(as.factor(dat5.agg$condition)), amount = 0.12)
dat5.agg$mean.resp2 <-
  jitter(dat5.agg$mean.resp, amount = 0.04)

boot.main_exp5_condition$ci.predicted$condition <-
  factor(boot.main_exp5_condition$ci.predicted$condition,
         levels = c("non-defeater", "defeater"))
ci_predicted_non_defeater <- boot.main_exp5_condition$ci.predicted %>%
  filter(condition == "non-defeater")
ci_predicted_defeater <- boot.main_exp5_condition$ci.predicted %>%
  filter(condition == "defeater")

# ggplot
exp5_plot_condition <-
  ggplot() +
  geom_point(
    data = dat5.agg,
    aes(x = condition2, y = mean.resp2, color = condition),
    size = 1.75, alpha = .4
  ) +
  scale_color_manual(values =
                       c("non-defeater" = "#4291F8", "defeater" = "#F19134")) +
  geom_line(
    data = dat5.agg,
    aes(x = condition2, y = mean.resp2, group = chimpanzee),
    color = "gray", lty = 1, alpha = .7
  ) +
  geom_errorbar(
    data = ci_predicted_non_defeater,
    aes(
      x = as.numeric(condition) - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = condition, color = condition
    ),
    width = 0.15, #size = 1
  ) +
  geom_errorbar(
    data = ci_predicted_defeater,
    aes(
      x = as.numeric(condition) - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      group = condition, color = condition
    ),
    width = 0.15, #size = 1
  ) +
  geom_point(
    data = ci_predicted_non_defeater,
    aes(x = as.numeric(condition) - 0.25, y = fitted),
    color = "#4291F8", size = 2
  ) +
  geom_point(
    data = ci_predicted_defeater,
    aes(x = as.numeric(condition) - 0.25, y = fitted),
    color = "#F19134", size = 2
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "gray"
  ) + # Adding the dashed line here
  scale_x_discrete(
    limits = c("non-defeater", "defeater"),
    name = "evidence condition",
    labels = c("non-defeater", "defeater")
  ) +
  scale_y_continuous(
    name = "proportion of belief revision",
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    limits = c(-0.05, 1),
    labels = scales::percent
  ) +
  geom_violin(
    data = dat5.agg,
    aes(x = condition, y = mean.resp, fill = condition),
    position = position_nudge(x = 0),
    alpha = .2
  ) +
  scale_fill_manual(values = c(
    "non-defeater" = "#4291F8",
    "defeater" = "#F19134"
  )) +
  labs(title = "Experiment 5") +
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

exp5_plot_condition

save.image("./R_images/study_5_analysis.RData")

