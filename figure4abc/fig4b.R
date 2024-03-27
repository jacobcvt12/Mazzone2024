library(tidyverse)
library(ggplot2)
library(patchwork)

# All fields included in this dataset are simulated for the purpose of code demonstration 
# and do not represent real results or actual data collected from participants.
data = readRDS("./data/cds_simulated_data.rds")

color_cv <- "#0072B2"
stage_color_pallete <- c("#8491B4FF", "#00A087FF", "#4DBBD5FF", "#E64B35FF")
font_size <- 5
n_iter <- 100000
expected_prop <- c(.546 / 0.999, .075 / 0.999, .218 / 0.999, .16 / 0.999)
cur_set = 'CV'

vars_list <- list(
  "HISTOLOGY" = c("Small Cell Carcinoma", "Non-Small Cell Carcinoma"),
  "TUMOR STAGE" = c(
    "T1", "T2", "T3", "T4",
    "T0, TX, Unk"
  ),
  "NODE STAGE" = c("N0", "N1", "N2", "N3", "NX, Unk"),
  "METASTASIS STAGE" = c("M0", "M1", "MX, Unk")
)

plot_data <- data[data$SET == cur_set, ]
sens_df <- tibble()
for (cur_var in names(vars_list)) {
  for (cur_row_name in vars_list[[cur_var]]) {
    stage_sens <- c()
    stages <- c()
    sens_low <- c()
    sens_up <- c()
    cur_data <- plot_data[plot_data[[cur_var]] == cur_row_name, ]
    
    if (nrow(cur_data) == 0) {
      next
    }
    
    # Calculate sensitivity per subgroup
    true_pos <- sum(cur_data[cur_data$CANCER_POS == TRUE, ]$DELFI_RESULT == "DLCST Positive")
    total_pos <- sum(cur_data$CANCER_POS)
    cutoff_sens <- true_pos / total_pos
    
    sens <- tibble(
      "Var" = cur_var,
      "Sub_var" = glue::glue("{cur_row_name} (n={nrow(cur_data)})"),
      "Cohort" = "CV",
      "Metric" = "Sensitivity",
      "Metric_val" = cutoff_sens,
      sens_low = binom::binom.wilson(
        x = true_pos,
        n = total_pos,
        conf.level = 0.95
      )[["lower"]],
      sens_up = binom::binom.wilson(
        x = true_pos,
        n = total_pos,
        conf.level = 0.95
      )[["upper"]]
    )
    
    sens_df <- dplyr::bind_rows(sens_df, sens)
  }
}

# Calculate Overall sensitivity
cur_data <- data[data$SET == cur_set, ]
true_pos <- sum(cur_data[cur_data$CANCER_POS == TRUE, ]$DELFI_RESULT == "DLCST Positive")
total_pos <- sum(cur_data$CANCER_POS)
cutoff_sens <- true_pos / total_pos

overall <- tibble(
  "Var" = "Overall",
  "Sub_var" = glue::glue("Observed (n={nrow(cur_data)})"),
  "Cohort" = "CV",
  "Metric" = "Sensitivity",
  "Metric_val" = cutoff_sens,
  sens_low = binom::binom.wilson(
    x = true_pos,
    n = total_pos,
    conf.level = 0.95
  )[["lower"]],
  sens_up = binom::binom.wilson(
    x = true_pos,
    n = total_pos,
    conf.level = 0.95
  )[["upper"]]
)

sens_df <- dplyr::bind_rows(overall, sens_df)
sens_df <- na.omit(sens_df)
sens_df$Var <- factor(sens_df$Var, levels = c(
  "HISTOLOGY",
  "TUMOR STAGE",
  "NODE STAGE",
  "METASTASIS STAGE",
  "STAGE",
  "Overall"
))
sens_df$Sub_var <- factor(
  sens_df$Sub_var,
  levels = sens_df$Sub_var
)

cohort_sens_df <- sens_df[sens_df$Metric == "Sensitivity", ]
cohort_sens_df <- cohort_sens_df[cohort_sens_df$Var != "Overall", ]

ggplot(
  data = cohort_sens_df,
  aes(
    x = Sub_var, y = Metric_val, ymin = sens_low, ymax = sens_up
  )
) +
  geom_hline(yintercept = overall$Metric_val, linetype = 1, size = .4) +
  geom_hline(yintercept = c(0, .2, .4, .6, .8, 1), linetype = 2, linewidth = .2) +
  geom_pointrange(fatten = .8, position = position_dodge(width = .5), colour = color_cv) +
  facet_grid(Metric ~ Var, scales = "free_x", space = "free") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 18), breaks = as.factor(cohort_sens_df$Sub_var)) +
  scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1)) +
  theme_classic() +
  ylab("Sensitivity") +
  xlab("") +
  theme(
    axis.text.y.right = element_blank(),
    axis.text.y = element_text(margin = margin(r = 2), size = font_size),
    panel.spacing = unit(2, "mm"),
    strip.background = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = .5),
    text = element_text(size = font_size),
    strip.text = element_text(size = font_size),
    strip.text.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = font_size)
  )
