############################################################################
# PLOTTING FIRST CHOICE
############################################################################
load("./R_images/additional_analyses.RData")
# Study Chimpanzee: Do Chimpanzees rationally revise their beliefs in light of
# evidence?
# Strength of evidence / First Choices in Study 1 and 2

# Last modified: Feb 27, 2024 12:50

############################################################################
# STUDY 1
############################################################################
library(tidyverse)
library(gghalves)
library(ggthemes)
library(cowplot)

xdata.agg <- xdata1 %>%
  group_by(Chimpanzee, first.evidence) %>%
  summarise(mean.resp = mean(first.evidence.chosen.numeric, na.rm = T)) %>%
  ungroup()
xdata.agg <- droplevels(xdata.agg)

xdata1$first.evidence <-
  factor(xdata1$first.evidence,
         levels = c(
           "visual_strong",
           "auditory_weak"
         )
  )
xdata.agg$first.evidence <-
  factor(xdata.agg$first.evidence,
         levels = c(
           "visual_strong",
           "auditory_weak"
         )
  )
levels(xdata.agg$first.evidence)

# Data manipulation outside ggplot
xdata.agg$first.evidence2 <-
  jitter(as.numeric(as.factor(xdata.agg$first.evidence)), amount = 0.12)
xdata.agg$mean.resp2 <-
  jitter(xdata.agg$mean.resp, amount = 0.04)

ci_predicted_study1_auditory <- boot.main_exp.1$ci.predicted %>%
  filter(first.evidence == "auditory_weak")
ci_predicted_study1_visual <- boot.main_exp.1$ci.predicted %>%
  filter(first.evidence == "visual_strong")

# ggplot
exp1_plot_first_choice <-
  ggplot() +
  geom_point(
    data = xdata.agg,
    aes(x = first.evidence2, y = mean.resp2, color = first.evidence),
    size = 2.5, alpha = .4
  ) +
  scale_color_manual(values = c(
    "visual_strong" = "dodgerblue",
    "auditory_weak" = "darkorange"
  )) +
  geom_line(
    data = xdata.agg,
    aes(x = first.evidence2, y = mean.resp2, group = Chimpanzee),
    color = "gray", lty = 1, alpha = .7
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray30") +
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
      "visual_strong",
      "auditory_weak"
    ),
    name = "type of evidence",
    labels = c(
      "visual (strong)",
      "auditroy (weak)"
    )
  ) +
  scale_y_continuous(
    name = "proportion of first choices",
    labels = scales::percent
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_violin(
    data = xdata.agg,
    aes(x = first.evidence, y = mean.resp, fill = first.evidence),
    position = position_nudge(x = 0),
    alpha = .2
  ) +
  scale_fill_manual(values = c(
    "visual_strong" = "dodgerblue",
    "auditory_weak" = "darkorange"
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
exp1_plot_first_choice

############################################################################
# PLOTTING STUDY 2
############################################################################
xdata.agg.2 <- xdata2 %>%
  mutate(first.evidence.chosen.numeric =
           as.numeric(as.factor(first.evidence.chosen)) - 1) %>%
  group_by(Chimpanzee, first.evidence) %>%
  summarise(mean.resp = mean(first.evidence.chosen.numeric, na.rm = T)) %>%
  ungroup()
xdata.agg.2 <- droplevels(xdata.agg.2)

xdata2$first.evidence <-
  factor(xdata2$first.evidence,
         levels = c(
           "auditory_strong",
           "trace_weak"
         )
  )
xdata.agg.2$first.evidence <-
  factor(xdata.agg.2$first.evidence,
         levels = c(
           "auditory_strong",
           "trace_weak"
         )
  )

# Data manipulation outside ggplot
xdata.agg.2$first.evidence2 <-
  jitter(as.numeric(as.factor(xdata.agg.2$first.evidence)), amount = 0.12)
xdata.agg.2$mean.resp2 <-
  jitter(xdata.agg.2$mean.resp, amount = 0.04)

ci_predicted_study2_trace <- boot.main_exp.2$ci.predicted %>%
  filter(first.evidence == "trace_weak")
ci_predicted_study2_auditory <- boot.main_exp.2$ci.predicted %>%
  filter(first.evidence == "auditory_strong")

# ggplot
exp2_plot_first_choice <-
  ggplot() +
  geom_point(
    data = xdata.agg.2,
    aes(x = first.evidence2, y = mean.resp2, color = first.evidence),
    size = 2.5, alpha = .4
  ) +
  scale_color_manual(values = c(
    "auditory_strong" = "dodgerblue",
    "trace_weak" = "darkorange"
  )) +
  geom_line(
    data = xdata.agg.2,
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
    data = ci_predicted_study2_auditory,
    aes(x = 1 - 0.25, y = fitted),
    color = "dodgerblue", size = 2.5
  ) +
  geom_point(
    data = ci_predicted_study2_trace,
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
    data = xdata.agg.2,
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

# Combine plots
library(ggpubr)
# theme_set(theme_pubr())
figure <- ggarrange(exp1_plot_first_choice, exp2_plot_first_choice,
                    labels = c(
                      "Experiment 1",
                      "Experiment 2"
                    ),
                    # label.x = 0,    # X position 0 for left
                    # label.y = 1,    # Y position 1 for top
                    hjust = -1, # hjust = 0 for left alignment
                    vjust = 2.5, # hjust = 0 for left alignment
                    ncol = 2, nrow = 1
)
figure
