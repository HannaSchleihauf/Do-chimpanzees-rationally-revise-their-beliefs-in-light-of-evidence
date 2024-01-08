# Study Name: Do chimpanzees rationally revise their beliefs in light of evidence?
# Stduy 1

# Last modified: Jan 8, 2024 21:17

############################################################################
# PACKAGES & FUNCTIONS
############################################################################
# library("readxl")
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

############################################################################
# DATA WRANGELING
############################################################################
xdata <-
  xdata %>%
  filter(!grepl('filler', Condition)) #deleted fillers

############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
xx.fe.re <- fe.re.tab(
  fe.model =
    "Belief_Revision ~ Condition*Trial",
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
  t.data$Trial - mean(t.data$Trial)

############################################################################
# FITTING THE MODEL AS PREREGISTERED
############################################################################
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000000))

full <- glmer(Belief_Revision ~
                Condition * Trial +
                (1 + Condition + Trial | Chimpanzee),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)


main <- glmer(Belief_Revision ~
                Condition + Trial +
                (1 + Condition + Trial | Chimpanzee),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

null <- glmer(Belief_Revision ~
                1 +
                (1 + Condition + Trial | Chimpanzee),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

summary(full)$varcor
ranef.diagn.plot(full)
round(summary(full)$coefficients, 3)

############################################################################
# CHECKING ASSUMPTIONS
############################################################################
overdisp.test(full)
library(car)
vif(main)
# Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full, contr = contr, use = c("Chimpanzee"))
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
round(anova(full, null, test = "Chisq"), 3)
round(drop1(full, test = "Chisq"), 3)
round(drop1(main, test = "Chisq"), 3)

plot(effect("Condition", main), type = "response")
plot(effect("Trial", main), type = "response")

# Coefficients of the full model
summary(full)$coefficients
