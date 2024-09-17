# Study Chimpanzee: Do Chimpanzees rationally revise their beliefs in light of
# evidence?
# Study 1

# Last modified: Feb 27, 2024 5:55

############################################################################
# PACKAGES & FUNCTIONS
############################################################################
library("easypackages")
libraries(
  "lme4", "tidyverse", "tidyselect", "optimx", "emmeans", "car",
  "effects", "agridat", "ggplot2", "ghibli", "ggdist", "see"
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
  ungroup() # Ungroup if needed for subsequent operations

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

tapply(
  as.numeric(xdata$Belief_Revision) - 1,
  list(xdata$Condition, as.factor(xdata$trial.per.Condition)), mean
)

ftable(xdata$Condition ~ xdata$Belief_Revision + xdata$trial.per.Condition)

############################################################################
# FITTING THE MODEL AS PREREGISTERED
############################################################################
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000000))

full_exp1 <- glmer(
  Belief_Revision ~
    Condition * z.Trial +
    (1 + Condition * z.Trial | Chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)

summary(full_exp1)$varcor
ranef.diagn.plot(full_exp1)
round(summary(full_exp1)$coefficients, 3) # massive SDs

# Model has a complete separation problem, thus we need to adjust the original
# data set

############################################################################
# MODEL SUFFERS FROM COMPLETE SEPERATION - ADJUSTED ANALYSIS
############################################################################
library("kyotil")

# check again which cells have a complete separation issue (have no 1s in it)
ftable(
  Belief_Revision ~ Condition + trial.per.Condition, t.data
) # Trial 2,3,4 in the strong_first condition
# determining in which cells we will be changing one data point at a time
to.change.1 <-
  (1:nrow(t.data))[t.data$Condition == "strong_first" &
    t.data$trial.per.Condition == "1"]
to.change.2 <-
  (1:nrow(t.data))[t.data$Condition == "strong_first" &
    t.data$trial.per.Condition == "2"]
to.change.3 <-
  (1:nrow(t.data))[t.data$Condition == "strong_first" &
    t.data$trial.per.Condition == "3"]
# load data storage object
load("./R_images/study_1_analysis_data_storage.RData")

# start with models
contr <-
  glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 10000000)
  )
for (i in 896:1000) { # i = 1
  set.seed(i)
  t.data$new.resp <- as.numeric(t.data$Belief_Revision) - 1
  xx <- sample(to.change.1, size = 1) # ample 1 element of to.change.1
  yy <- sample(to.change.2, size = 1)
  zz <- sample(to.change.3, size = 1)
  all.res$to.change.1[i] <- xx # row number into data storage
  all.res$to.change.2[i] <- yy
  all.res$to.change.3[i] <- zz
  t.data$new.resp[xx] <- ifelse(t.data$new.resp[xx] == 0, 1, 0)
  t.data$new.resp[yy] <- ifelse(t.data$new.resp[yy] == 0, 1, 0)
  t.data$new.resp[zz] <- ifelse(t.data$new.resp[zz] == 0, 1, 0)
  # Full model
  full1 <-
    keepWarnings(
      glmer(
        new.resp ~
          Condition * z.Trial +
          (1 + Condition * z.Trial | Chimpanzee),
        data = t.data, control = contr,
        family = binomial(link = "logit")
      )
    )

  red1 <-
    keepWarnings(
      glmer(
        new.resp ~
          Condition + z.Trial +
          (1 + Condition * z.Trial | Chimpanzee),
        data = t.data, control = contr,
        family = binomial(link = "logit")
      )
    )

  null1 <-
    keepWarnings(
      glmer(
        new.resp ~
          1 +
          (1 + Condition * z.Trial | Chimpanzee),
        data = t.data, control = contr,
        family = binomial(link = "logit")
      )
    )

  if (length(full1$warnings) == 0 &
    length(red1$warnings) == 0 &
    length(null1$warnings) == 0) {
    full <- full1$value
    red <- red1$value
    null <- null1$value

    all.res$warnings[i] <- "no"

    all.res$all.full.coeffs.intercept[i] <-
      summary(full)$coefficients["(Intercept)", 1]
    all.res$all.full.coeffs.Condition[i] <-
      summary(full)$coefficients["Conditionweak_first", 1]
    all.res$all.full.coeffs.Trial[i] <-
      summary(full)$coefficients["z.Trial", 1]
    all.res$all.full.coeffs.Condition.Trial[i] <-
      summary(full)$coefficients["Conditionweak_first:z.Trial", 1]

    # full null model comparisons
    test.full.null <-
      as.data.frame(anova(null, full, test = "Chisq"))[
        "full",
        c("Chisq", "Df", "Pr(>Chisq)")
      ]
    all.res$test.full.null.Chisq[i] <- test.full.null$Chisq
    all.res$test.full.null.Df[i] <- test.full.null$Df
    all.res$test.full.null.Pr..Chisq[i] <- test.full.null$`Pr(>Chisq)`

    # full red model comparisons
    tests.full.red <-
      drop1p(
        model.res = full, para = F, contr = contr,
        n.cores = c("all-1", "all"), to.del = NULL
      )
    test.2.way.int <-
      as.data.frame(tests.full.red$drop1.res[2, c(
        "Chisq", "Chi.Df", "Pr..Chisq.", "n.opt.warnings", "n.fun.warnings"
      )])

    all.res$test.2.way.Condition.Trial.Chisq[i] <-
      test.2.way.int$Chisq
    all.res$test.2.way.Condition.Trial.Chi.Df[i] <-
      test.2.way.int$Chi.Df
    all.res$test.2.way.Condition.Trial.Pr..Chisq[i] <-
      test.2.way.int$Pr..Chisq.
    all.res$test.2.way.Condition.Trial.n.opt.warnings[i] <-
      test.2.way.int$n.opt.warnings
    all.res$test.2.way.Condition.Trial.n.fun.warnings[i] <-
      test.2.way.int$n.fun.warnings

    # assumptions
    # colliniarity
    all.res$colliniarity.test.Condition[i] <-
      vif(red)[1]
    all.res$colliniarity.test.Trial[i] <-
      vif(red)[2]

    # red-main model comparisons
    tests.red.main <-
      drop1p(
        model.res = red, para = F, contr = contr,
        n.cores = c("all-1", "all"), to.del = NULL
      )
    test.main.Condition <-
      as.data.frame(tests.red.main$drop1.res[2, c(
        "Chisq", "Chi.Df", "Pr..Chisq.", "n.opt.warnings", "n.fun.warnings"
      )])
    test.main.Trial <-
      as.data.frame(tests.red.main$drop1.res[3, c(
        "Chisq", "Chi.Df", "Pr..Chisq.", "n.opt.warnings", "n.fun.warnings"
      )])
    all.res$test.main.Condition.Chisq[i] <-
      test.main.Condition$Chisq
    all.res$test.main.Condition.Chi.Df[i] <-
      test.main.Condition$Chi.Df
    all.res$test.main.Condition.Pr..Chisq[i] <-
      test.main.Condition$Pr..Chisq.
    all.res$test.main.Condition.n.opt.warnings[i] <-
      test.main.Condition$n.opt.warnings
    all.res$test.main.Condition.n.fun.warnings[i] <-
      test.main.Condition$n.fun.warnings
    all.res$test.main.Trial.Chisq[i] <-
      test.main.Trial$Chisq
    all.res$test.main.Trial.Chi.Df[i] <-
      test.main.Trial$Chi.Df
    all.res$test.main.Trial.Pr..Chisq[i] <-
      test.main.Trial$Pr..Chisq.
    all.res$test.main.Trial.n.opt.warnings[i] <-
      test.main.Trial$n.opt.warnings
    all.res$test.main.Trial.n.fun.warnings[i] <-
      test.main.Trial$n.fun.warnings

    # boot full model
    boot.full <- boot.glmm.pred(
      model.res = full, excl.warnings = T, nboots = 10,
      para = F, resol = 100, level = 0.95, use = c("Condition", "z.Trial")
    )
    boot.full$ci.predicted$name <-
      paste(boot.full$ci.predicted$Condition,
        boot.full$ci.predicted$Z.Trial,
        sep = "."
      )

    boot.full.values <-
      rbind(boot.full.values, data.frame(
        term = boot.full$ci.predicted$Condition,
        fitted = boot.full$ci.predicted$fitted,
        lower.cl = boot.full$ci.predicted$lower.cl,
        upper.cl = boot.full$ci.predicted$upper.cl
      ))

    boot.full.estimates <-
      rbind(boot.full.estimates, data.frame(
        term = rownames(boot.full$ci.estimates),
        orig = boot.full$ci.estimates$orig,
        X2.5. = boot.full$ci.estimates$X2.5.,
        X97.5. = boot.full$ci.estimates$X97.5.
      ))

    # boot Cond model
    boot.Cond <- boot.glmm.pred(
      model.res = full, excl.warnings = T, nboots = 10,
      para = F, resol = 100, level = 0.95, use = c("Condition")
    )
    boot.Cond.values <-
      rbind(boot.Cond.values, data.frame(
        term = boot.Cond$ci.predicted$Condition,
        fitted = boot.Cond$ci.predicted$fitted,
        lower.cl = boot.Cond$ci.predicted$lower.cl,
        upper.cl = boot.Cond$ci.predicted$upper.cl
      ))
    boot.Cond.estimates <-
      rbind(boot.Cond.estimates, data.frame(
        term = rownames(boot.Cond$ci.estimates),
        orig = boot.Cond$ci.estimates$orig,
        X2.5. = boot.Cond$ci.estimates$X2.5.,
        X97.5. = boot.Cond$ci.estimates$X97.5.
      ))
  } else {
    all.res$warnings[i] <- "yes"
  }
}

save.image("./R_images/study_1_analysis_complete_separation.RData")

############################################################################
# EVALUATION OF THE RESULTS
############################################################################

# how many models did converge
sum(all.res$warnings == "no")

# means of full-null-comparisons
# Chisq
round(mean(all.res$test.full.null.Chisq, na.rm = T), 10)
round(range(all.res$test.full.null.Chisq, na.rm = T), 10)
# DF
range(all.res$test.full.null.Df, na.rm = T)
# p-value
round(mean(all.res$test.full.null.Pr..Chisq, na.rm = T), 10)
round(range(all.res$test.full.null.Pr..Chisq, na.rm = T), 10)

# colliniarity
round(mean(all.res$colliniarity.test.Condition, na.rm = T), 10)
round(range(all.res$colliniarity.test.Condition, na.rm = T), 10)
round(mean(all.res$colliniarity.test.Trial, na.rm = T), 10)
round(range(all.res$colliniarity.test.Trial, na.rm = T), 10)

# means of two-way interaction
# Chisq
round(mean(all.res$test.2.way.Condition.Trial.Chisq, na.rm = T), 10)
round(range(all.res$test.2.way.Condition.Trial.Chisq, na.rm = T), 10)
# DF
range(all.res$test.2.way.Condition.Trial.Chi.Df, na.rm = T)
# p-value
round(mean(all.res$test.2.way.Condition.Trial.Pr..Chisq, na.rm = T), 10)
round(range(all.res$test.2.way.Condition.Trial.Pr..Chisq, na.rm = T), 10)

# Main Effect Trial
# Chisq
round(mean(all.res$test.main.Trial.Chisq, na.rm = T), 10)
round(range(all.res$test.main.Trial.Chisq, na.rm = T), 10)
# DF
range(all.res$test.main.Trial.Chi.Df, na.rm = T)
# p-value
round(mean(all.res$test.main.Trial.Pr..Chisq, na.rm = T), 10)
round(range(all.res$test.main.Trial.Pr..Chisq, na.rm = T), 10)

# Main Effect Condition
# Chisq
round(mean(all.res$test.main.Condition.Chisq, na.rm = T), 10)
round(range(all.res$test.main.Condition.Chisq, na.rm = T), 10)
# DF
range(all.res$test.main.Condition.Chi.Df, na.rm = T)
# p-value
round(mean(
  all.res$test.main.Condition.Pr..Chisq,
  na.rm = T
), 10)
round(range(all.res$test.main.Condition.Pr..Chisq, na.rm = T), 10)

# coefs of the full  model
all.coefs <-
  all.res %>%
  select(vars_select(
    names(all.res),
    starts_with("all.full", ignore.case = TRUE)
  ))

round(colMeans(all.coefs, na.rm = T), 3)
data.frame(
  min = sapply(all.coefs, min, na.rm = T),
  max = sapply(all.coefs, max, na.rm = T)
)

# boot.full fitted values and confidence intervals
boot.full <- mapply(
  FUN = tapply, X = boot.full.values[, -1],
  MoreArgs = list(INDEX = boot.full.values$term, FUN = mean)
)
results.boot.full.values <- round(boot.full, 3)

# boot.full estimates
xx <- mapply(
  FUN = tapply, X = boot.full.estimates[, -1],
  MoreArgs = list(INDEX = boot.full.estimates$term, FUN = mean)
)
results.boot.full.estimates <- round(xx, 3)

# boot.Cond fitted values and confidence intervals
boot.Cond <- mapply(
  FUN = tapply, X = boot.Cond.values[, -1],
  MoreArgs = list(INDEX = boot.Cond.values$term, FUN = mean)
)
results.boot.Cond.values <- round(boot.Cond, 3)

# boot.Cond estimates
xx <- mapply(
  FUN = tapply, X = boot.Cond.estimates[, -1],
  MoreArgs = list(INDEX = boot.Cond.estimates$term, FUN = mean)
) # median
results.boot.Cond.estimates <- round(xx, 3)




# Code to test the random effects of 10 randomly selected models
i <- 978

xx <- all.res$to.change.1[i]
yy <- all.res$to.change.2[i]
zz <- all.res$to.change.3[i]
t.data$new.resp[xx] <- ifelse(t.data$new.resp[xx] == 0, 1, 0)
t.data$new.resp[yy] <- ifelse(t.data$new.resp[yy] == 0, 1, 0)
t.data$new.resp[zz] <- ifelse(t.data$new.resp[zz] == 0, 1, 0)

full_exp1 <- glmer(
  new.resp ~
    Condition * z.Trial +
    (1 + Condition * z.Trial | Chimpanzee),
  data = t.data, control = contr,
  family = binomial(link = "logit")
)

summary(full_exp1)$varcor
ranef.diagn.plot(full_exp1)
round(summary(full_exp1)$coefficients, 3)

############################################################################
# PLOTTING
############################################################################
load("./R_images/study_1_analysis_complete_separation.RData")

xdata.agg <- xdata %>%
  mutate(belief_revision.numeric = as.numeric(Belief_Revision) - 1) %>%
  group_by(Chimpanzee, Condition) %>%
  summarise(mean.resp = mean(belief_revision.numeric, na.rm = T)) %>%
  ungroup()
xdata.agg$Condition <- droplevels(xdata.agg$Condition)

xdata.agg <- xdata.agg %>%
  group_by(Condition, mean.resp) %>%
  mutate(trial.per.cell = row_number()) %>%
  ungroup() %>%
  mutate(Condition2 = as.numeric(as.factor(Condition)) -
    0.08 - (trial.per.cell * 0.04))

results.boot.Cond.values <- data.frame(results.boot.Cond.values)
results.boot.Cond.values$Condition <- rownames(results.boot.Cond.values)
rownames(results.boot.Cond.values) <- NULL

# ggplot
exp1_plot_Condition <-
  ggplot() +
  geom_line(
    data = xdata.agg,
    aes(
      x = as.numeric(Condition2),
      y = mean.resp,
      group = Chimpanzee
    ),
    color = "gray", lty = 1, alpha = 0.7
  ) +
  geom_errorbar(
    data = results.boot.Cond.values,
    aes(
      x = as.numeric(Condition),
      y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      color = Condition # Color should now apply properly
    ),
    width = 0.1, size = 1, color = c("#4291F8", "#F19134")
  ) +
  geom_point(
    data = results.boot.Cond.values,
    aes(
      x = as.numeric(Condition), y = fitted,
      color = Condition
    ), size = 2.5, color = c("#4291F8", "#F19134")
  ) +
  geom_point(
    data = xdata.agg,
    aes(x = Condition2, y = mean.resp, color = Condition),
    size = 2.5, alpha = .4
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
  geom_violinhalf(
    data = xdata.agg,
    aes(x = Condition, y = mean.resp, fill = Condition),
    position = position_nudge(x = 0.1),
    alpha = .2
  ) +
  scale_fill_manual(values = c(
    "strong_first" = "#4291F8",
    "weak_first" = "#F19134"
  )) +
  scale_color_manual(
    values =
      c("strong_first" = "#4291F8", "weak_first" = "#F19134")
  ) +
  labs(title = "Experiment 1") +
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
exp1_plot_Condition

save.image("./R_images/study_1_analysis_complete_separation.RData")
