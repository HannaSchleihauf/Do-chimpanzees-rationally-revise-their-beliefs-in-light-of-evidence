##################### Prepare data for model fitting ###########################
secondChoiceData_exp1_GLM <- dat1 %>%
  filter(!grepl("filler", Condition))  %>% # deleted fillers
  arrange(Chimpanzee, Trial) %>%
  group_by(Condition) %>%
  mutate(trial.per.Condition = row_number()) %>%
  mutate(second.choice.numeric = case_when(
    Condition == "strong_first" & Second_Choice == "strong" ~ "1",
    Condition == "strong_first" & Second_Choice == "weak" ~ "0",
    Condition == "weak_first" & Second_Choice == "weak" ~ "0",
    Condition == "weak_first" & Second_Choice == "strong" ~ "1"
  ))

secondChoiceData_exp1_GLM$first.evidence.chosen.numeric <-
  as.numeric(secondChoiceData_exp1_GLM$first.evidence.chosen.numeric)

# Table
second.exp1.fe.re <- fe.re.tab(fe.model =
                                 "second.choice.numeric ~ Condition + trial.per.Condition",
                               re = "(1|Chimpanzee)", data = secondChoiceData_exp1_GLM)
second.exp1.t.data <- second.exp1.fe.re$data
str(second.exp1.t.data)

# Center/scale dummy variables
second.exp1.t.data$Condition.weak_first.code <-
  second.exp1.t.data$Condition.weak_first -
  mean(second.exp1.t.data$Condition.weak_first)
second.exp1.t.data$z.Trial <- scale(second.exp1.t.data$trial.per.Condition)


######################### Fit preregistered GLM ################################

# define controls
contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000000))

# fit three model levels
second_full_exp.1 <- glmer(second.choice.numeric ~ Condition * z.Trial +
                             (1 + Condition.weak_first.code *
                                z.Trial || Chimpanzee),
                           data = second.exp1.t.data, control = contr,
                           family = binomial(link = "logit")
)

second_main_exp.1 <- glmer(second.choice.numeric ~ Condition + z.Trial +
                             (1 + Condition.weak_first.code *
                                z.Trial || Chimpanzee),
                           data = second.exp1.t.data, control = contr,
                           family = binomial(link = "logit")
)

second_null_exp.1 <- glmer(second.choice.numeric ~ 1 +
                             (1 + Condition.weak_first.code *
                                z.Trial || Chimpanzee),
                           data = second.exp1.t.data, control = contr,
                           family = binomial(link = "logit")
)

summary(second_null_exp.1) # This shows that it is different from chance!

# Model diagnostics
summary(second_full_exp.1)$varcor
ranef.diagn.plot(second_full_exp.1)
round(summary(second_full_exp.1)$coefficients, 3)


######################## Check model assumptions ###############################

# overdispersion
overdisp.test(second_full_exp.1)

# variance inflation factor
vif(second_main_exp.1)

# model stability
m.stab.b <- glmm.model.stab(model.res = second_full_exp.1,
                            contr = contr, use = c("Chimpanzee"))
m.stab.b$detailed$warnings
as.data.frame(round(m.stab.b$summary[, -1], 3))
m.stab.plot(round(m.stab.b$summary[, -1], 3))

# model stability only with models that did not give a warning message
m.stab.b$detailed <- m.stab.b$detailed %>%
  filter(warnings == "none")
cbind(apply(m.stab.b$detailed, 2, min),apply(m.stab.b$detailed, 2, max))

########################### Model Comparisons ##################################

round(anova(second_full_exp.1, second_null_exp.1, test = "Chisq"), 3)
round(drop1(second_full_exp.1, test = "Chisq"), 3)
round(drop1(second_main_exp.1, test = "Chisq"), 3)

# Coefficients of the full_exp.1 model
round(summary(second_full_exp.1)$coefficients, 3)
library("effects")
## First peek at effects
plot(effect("Condition", second_main_exp.1), type = "response")
plot(effect("z.Trial", main_exp1), type = "response")


############################# Bootstraps #######################################
# already run and saved in image

# For second_full_exp.1 model
boot.second_full_exp.1 <- boot.glmm.pred(model.res = second_full_exp.1,
                                         excl.warnings = T,
                                         nboots = 1000, para = F, level = 0.95,
                                         use = c("second.evidence", "z.Trial"))
# Save predicted CIs
exp1_second_full_ci_predicted <- boot.second_full_exp.1$ci.predicted

# Plot
m.stab.plot(round(boot.second_full_exp.1$ci.estimates, 3))

# For main_exp.1 model
boot.second_main_exp.1 <- boot.glmm.pred(model.res = second_main_exp.1,
                                         excl.warnings = T,
                                         nboots = 1000, para = F, level = 0.95,
                                         use = c("second.evidence"))

# Save predicted CIs
exp1_second_main_ci_predicted <- boot.second_main_exp.1$ci.predicted

# Plot
m.stab.plot(round(boot.second_main_exp.1$ci.estimates, 3))
```
