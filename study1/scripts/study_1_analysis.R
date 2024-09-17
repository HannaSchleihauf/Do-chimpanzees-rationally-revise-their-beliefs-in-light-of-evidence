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
  group_by(Chimpanzee, Condition) %>%
  mutate(trial.per.Condition = row_number()) %>%
  ungroup()  # Ungroup if needed for subsequent operations

str(xdata)
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

tapply(as.numeric(xdata$Belief_Revision)-1,
       list(xdata$Condition, as.factor(xdata$trial.per.Condition)), mean)

ftable(xdata$Condition ~ xdata$Belief_Revision + xdata$trial.per.Condition)

############################################################################
# FITTING THE MODEL AS PREREGISTERED
############################################################################
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000000))

full_exp1 <- glmer(Belief_Revision ~
  Condition * z.Trial +
  (1 + Condition * z.Trial | Chimpanzee),
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

# Model has a complete separation problem.
