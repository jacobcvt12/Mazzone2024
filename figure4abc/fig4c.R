library(tidyverse)
library(ggplot2)
library(patchwork)

# Data is simulated for the purpose of code demonstration. It doesn't represent the real data.
data = readRDS("./data/cds_simulated_data.rds")

bootstrap_ci <- function(df, iter = 10000) {
  sens_list <- c()
  size <- sum(df$CANCER_POS)
  x <- sum(df[df$CANCER_POS == TRUE, ]$DELFI_RESULT == "DLCST Positive")
  for (i in 1:iter) {
    cur_sample <- rbinom(n = size, size = 1, prob = x / size)
    cur_sens <- mean(cur_sample)
    sens_list <- c(sens_list, cur_sens)
  }
  return(sens_list)
}

color_cv <- "#0072B2"
stage_color_pallete <- c("#8491B4FF", "#00A087FF", "#4DBBD5FF", "#E64B35FF")
font_size <- 5
n_iter <- 100000
expected_prop <- c(.546 / 0.999, .075 / 0.999, .218 / 0.999, .16 / 0.999)
cur_set = 'CV'

set.seed(123)

# Calculate sensitivities by stage
sens_df <- tibble()
stage_sens <- c()
stages <- c()
sens_low <- c()
sens_up <- c()
for (cur_stage in c("I", "II", "III", "IV")) {
  cur_data <- data[data$SET == cur_set & data$STAGE == cur_stage, ]

  if (nrow(data) == 0) {
    next
  }

  cur_df <- cur_data[cur_data$LCNSTAGE %in% c(cur_stage), ]

  if (nrow(cur_df) == 0) {
    cur_cutoff_sens <- NA
    cur_sens_low <- 0
    cur_sens_up <- 0
  } else {
    cur_cutoff_sens <- sum(cur_df[cur_df$CANCER_POS == TRUE, ]$DELFI_RESULT == "DLCST Positive") / (sum(cur_df$CANCER_POS))

    cur_sens_low <- binom::binom.wilson(
      x = sum(cur_df[cur_df$CANCER_POS == TRUE, ]$DELFI_RESULT == "DLCST Positive"),
      n = sum(cur_df$CNCLBL == "Lung Cancer"),
      conf.level = 0.95
    )[["lower"]]

    cur_sens_up <- binom::binom.wilson(
      x = sum(cur_df[cur_df$CANCER_POS == TRUE, ]$DELFI_RESULT == "DLCST Positive"),
      n = sum(cur_df$CNCLBL == "Lung Cancer"),
      conf.level = 0.95
    )[["upper"]]
  }

  sens <- tibble(
    "Var" = "Stage",
    "Cat" = cur_stage,
    "Sub_var" = glue::glue("{cur_stage} (n={nrow(cur_df)})"),
    "Cohort" = cur_set,
    "Metric" = "Sensitivity",
    "Metric_val" = cur_cutoff_sens,
    sens_low = cur_sens_low,
    sens_up = cur_sens_up
  )

  sens_df <- dplyr::bind_rows(sens_df, sens)
}

# Calculate overall sensitivity
cur_data <- data[data$SET == cur_set, ]
cur_cutoff_sens <- sum(cur_data[cur_data$CANCER_POS == TRUE, ]$DELFI_RESULT == "DLCST Positive") / (sum(cur_data$CANCER_POS))

cur_sens_low <- binom::binom.wilson(
  x = sum(cur_data[cur_data$CANCER_POS == TRUE, ]$DELFI_RESULT == "DLCST Positive"),
  n = sum(cur_data$CNCLBL == "Lung Cancer"),
  conf.level = 0.95
)[["lower"]]

cur_sens_up <- binom::binom.wilson(
  x = sum(cur_data[cur_data$CANCER_POS == TRUE, ]$DELFI_RESULT == "DLCST Positive"),
  n = sum(cur_data$CNCLBL == "Lung Cancer"),
  conf.level = 0.95
)[["upper"]]

sens <- tibble(
  "Var" = "Overall",
  "Cat" = "Overall",
  "Sub_var" = glue::glue("Overall (n={nrow(cur_data)})"),
  "Cohort" = cur_set,
  "Metric" = "Sensitivity",
  "Metric_val" = cur_cutoff_sens,
  sens_low = cur_sens_low,
  sens_up = cur_sens_up
)

sens_df <- dplyr::bind_rows(sens, sens_df)

# Calculate overall, stage-weighted sensitivity
stage_sens <- c()
stages <- c()
sens_low <- c()
sens_up <- c()
stage_sens <- matrix(data = NA, nrow = n_iter, ncol = 4)
i <- 1
for (cur_stage in c("I", "II", "III", "IV")) {
  stages <- c(stages, cur_stage)
  if (nrow(cur_data) == 0) {
    next
  }

  cur_df <- cur_data[cur_data$LCNSTAGE %in% c(cur_stage), ]
  stg_sens_list <- bootstrap_ci(cur_df, iter = n_iter)
  stage_sens[, i] <- stg_sens_list
  i <- i + 1

  if (nrow(cur_df) == 0) {
    next
  }
}

stg_wtd_mean_overall <- stage_sens %*% expected_prop
sens_mean <- mean(stg_wtd_mean_overall)
stg_weighted_sens_ci <- quantile(stg_wtd_mean_overall, probs = c(.025, .975))
sens_low <- stg_weighted_sens_ci[1]
sens_up <- stg_weighted_sens_ci[2]
stg_wtd_sens_overall <- tibble(
  "Var" = "Overall",
  "Cat" = "Stage Weighted",
  "Sub_var" = glue::glue("Stage Weighted (n={nrow(cur_data)})"),
  "Cohort" = cur_set,
  "Metric" = "Sensitivity",
  "Metric_val" = sens_mean,
  "sens_low" = sens_low,
  "sens_up" = sens_up
)

sens_df <- dplyr::bind_rows(sens_df, stg_wtd_sens_overall)

sens_df <- na.omit(sens_df)

sens_df$Sub_var <- factor(sens_df$Sub_var, levels = sens_df$Sub_var)

# Plot
cohort_sens_df <- sens_df[sens_df$Cohort == cur_set, ]
stage_cohort_sens_df <- cohort_sens_df[!cohort_sens_df$Cat %in% c("Overall", "Stage Weighted"), ]
stage_cohort_sens_df$Props <- c(expected_prop)

g1 <- ggplot(
  data = stage_cohort_sens_df,
  aes(
    x = Sub_var, y = Metric_val, ymin = sens_low, ymax = sens_up,
    color = Sub_var
  )
) +
  geom_hline(yintercept = c(0, .2, .4, .6, .8, 1), linetype = 2, linewidth = .2) +
  geom_pointrange(fatten = 3, position = position_dodge(width = .5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 18), breaks = as.factor(cohort_sens_df$Sub_var)) +
  scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1)) +
  ylab("Sensitivity") +
  xlab("") +
  ggtitle("Sensitivity") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = font_size),
    axis.text.y = element_text(margin = margin(r = 2), size = font_size),
    panel.spacing = unit(2, "mm"),
    strip.background = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = .2),
    text = element_text(size = font_size),
    strip.text = element_text(size = font_size),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = font_size)
  ) +
  scale_color_manual(values = stage_color_pallete)

stage_cohort_sens_df <- stage_cohort_sens_df[stage_cohort_sens_df$Var != "Overall", ]
stage_cohort_sens_df$Label <- c(
  "I, 54.7%", "II, 7.5%",
  "III, 21.8%", "IV, 16.0%"
)

stage_cohort_sens_df$Sub_var <- factor(stage_cohort_sens_df$Sub_var, levels = c(
  "IV",
  "III",
  "II",
  "I"
))
stage_cohort_sens_df$Cat <- factor(stage_cohort_sens_df$Cat, levels = c(
  "IV",
  "III",
  "II",
  "I"
))

g2 <- ggplot(
  data = stage_cohort_sens_df,
  aes(
    x = Var, y = Props, fill = Cat, label = Label
  )
) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(size = 1.75, position = position_stack(vjust = 0.5)) +
  scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1)) +
  ylab("Prevalence") +
  xlab("") +
  ggtitle("Prevalence in \nscreening population") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = font_size),
    axis.text.y = element_text(margin = margin(r = 2), size = font_size),
    axis.text.x = element_text(size = font_size),
    panel.spacing = unit(2, "mm"),
    strip.background = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = .2),
    text = element_text(size = font_size),
    strip.text = element_text(size = font_size),
    legend.position = "none"
  ) +
  scale_fill_manual(values = rev(stage_color_pallete))

g3 <- ggplot(
  data = cohort_sens_df[cohort_sens_df$Cat == "Stage Weighted", ],
  aes(
    x = Sub_var, y = Metric_val, ymin = sens_low, ymax = sens_up,
    color = "black"
  )
) +
  geom_hline(yintercept = cohort_sens_df[cohort_sens_df$Cat == "Stage Weighted", ]$Metric_val, linetype = 1, size = .4) +
  geom_hline(yintercept = c(0, .2, .4, .6, .8, 1), linetype = 2, linewidth = .2) +
  geom_pointrange(fatten = 3, position = position_dodge(width = .5), color = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 18), breaks = as.factor(cohort_sens_df$Sub_var)) +
  scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1)) +
  ylab("Sensitivity") +
  xlab("") +
  ggtitle("Sensitivity in \nscreening population") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = font_size),
    axis.text.y = element_text(margin = margin(r = 2), size = font_size),
    axis.text.x = element_text(size = font_size),
    panel.spacing = unit(2, "mm"),
    strip.background = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = .2),
    text = element_text(size = font_size),
    strip.text = element_text(size = font_size),
    legend.position = "none"
  )

g1 + g2 + plot_spacer() + g3 +
  plot_layout(widths = c(3, 1, .5, 1))
