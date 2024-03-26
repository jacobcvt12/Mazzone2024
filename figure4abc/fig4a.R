library(binom)
library(dplyr)
library(ggplot2)
library(glue)
library(purrr)
library(rlang)
library(stringr)
library(tibble)

# Data is simulated for the purpose of code demonstration. It doesn't represent the real data.
data = readRDS("./data/cds_simulated_data.rds")

point_size <- 0.25
font_size_x <- 5
font_size_y <- 5

cbCol <- c("#0072B2", "#E69F00", "#009E73", "#F0E473", "#56B4E9", "#D55E00", "#CC79A7", "#333333")


calculate_accuracy <- function(data, cat_var, type = c("sens", "spec")) {
  labels <- switch(type,
                   "sens" = c("DLCST Positive", "Lung Cancer"),
                   "spec" = c("DLCST Negative", "Non-cancer")
  )
  accuracy <- data %>%
    dplyr::group_by(!!rlang::sym(cat_var)) %>%
    dplyr::summarise(
      correct_call = sum(DELFI_RESULT == labels[[1]] & CNCLBL == labels[[2]]),
      cases = sum(CNCLBL == labels[[2]]),
      total = n(),
      accuracy = correct_call / cases,
      ci_lower = binom::binom.wilson(x = correct_call, n = cases, conf.level = 0.95)[["lower"]],
      ci_upper = binom::binom.wilson(x = correct_call, n = cases, conf.level = 0.95)[["upper"]]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(cases, correct_call)) %>%
    dplyr::mutate(label = type)
  return(accuracy)
}

combine_data <- function(cat_var, var_lab, ...) {
  data <- list(...) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    dplyr::group_by(!!rlang::sym(cat_var)) %>%
    dplyr::mutate(
      order = factor(dplyr::row_number()),
      clinical_var = var_lab,
      label = factor(label,
                     levels = c("sens", "spec"),
                     labels = c("Sensitivity", "Specificity")
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(category = !!rlang::sym(cat_var))
  return(data)
}

combine_vars <- function(data) {
  overall_sens <- calculate_accuracy(data, "", type = "sens")
  overall_spec <- calculate_accuracy(data, "", type = "spec")
  overall <- combine_data(
    cat_var = "", var_lab = "Overall",
    overall_sens, overall_spec
  ) %>%
    dplyr::mutate(category = "Overall")
  
  age_cat_sens <- calculate_accuracy(data, "AGECAT", type = "sens")
  age_cat_spec <- calculate_accuracy(data, "AGECAT", type = "spec")
  age_cat <- combine_data(
    cat_var = "AGECAT", var_lab = glue::glue("Age\nBinary"),
    age_cat_sens, age_cat_spec
  )
  
  sex_sens <- calculate_accuracy(data, "SEX", "sens")
  sex_spec <- calculate_accuracy(data, "SEX", "spec")
  sex <- combine_data(
    cat_var = "SEX", var_lab = "Sex",
    sex_sens, sex_spec
  )
  
  race_new_sens <- calculate_accuracy(data, "RACE_NEW", "sens")
  race_new_spec <- calculate_accuracy(data, "RACE_NEW", "spec")
  race_new <- combine_data(
    cat_var = "RACE_NEW", var_lab = "Race",
    race_new_sens, race_new_spec
  )
  
  eth_sens <- calculate_accuracy(data, "ETHNIC_NEW", "sens")
  eth_spec <- calculate_accuracy(data, "ETHNIC_NEW", "spec")
  eth <- combine_data(
    cat_var = "ETHNIC_NEW", var_lab = "Ethnicity",
    eth_sens, eth_spec
  )
  
  data_copd <- data %>%
    dplyr::filter(MH_COPD == "COPD")
  copd_sens <- calculate_accuracy(data_copd, "MH_COPD", "sens")
  copd_spec <- calculate_accuracy(data_copd, "MH_COPD", "spec")
  copd <- combine_data(
    cat_var = "MH_COPD", var_lab = "Comorbidities",
    copd_sens, copd_spec
  )
  
  data_diab <- data %>%
    dplyr::filter(MH_DIAB == "Type-II Diabetes")
  diab_sens <- calculate_accuracy(data_diab, "MH_DIAB", "sens")
  diab_spec <- calculate_accuracy(data_diab, "MH_DIAB", "spec")
  diab <- combine_data(
    cat_var = "MH_DIAB", var_lab = "Comorbidities",
    diab_sens, diab_spec
  )
  
  data_bmi <- data %>%
    dplyr::filter(BMI_CAT == "Overweight/Obesity, BMI > 25")
  bmi_sens <- calculate_accuracy(data_bmi, "BMI_CAT", "sens")
  bmi_spec <- calculate_accuracy(data_bmi, "BMI_CAT", "spec")
  bmi <- combine_data(
    cat_var = "BMI_CAT", var_lab = "Comorbidities",
    bmi_sens, bmi_spec
  )

  pch_sens <- calculate_accuracy(data, "PRCNRYN", "sens")
  pch_spec <- calculate_accuracy(data, "PRCNRYN", "spec")
  pch <- combine_data(
    cat_var = "PRCNRYN", var_lab = glue::glue("Prior\nCancer"),
    pch_sens, pch_spec
  )
  
  smoking_sens <- calculate_accuracy(data, "SUCAT", "sens")
  smoking_spec <- calculate_accuracy(data, "SUCAT", "spec")
  smoking <- combine_data(
    cat_var = "SUCAT", var_lab = glue::glue("Smoking\nStatus"),
    smoking_sens, smoking_spec
  )
  
  packyrs_sens <- calculate_accuracy(data, "SUPACKYR_CAT", "sens")
  packyrs_spec <- calculate_accuracy(data, "SUPACKYR_CAT", "spec")
  packyrs <- combine_data(
    cat_var = "SUPACKYR_CAT", var_lab = glue::glue("Pack-\nYears"),
    packyrs_sens, packyrs_spec
  )
  
  combined <- list(
    overall, age_cat, sex, race_new,
    eth, copd, diab, bmi, pch, 
    smoking, packyrs
  ) %>%
    purrr::reduce(dplyr::bind_rows)
  
  return(combined)
}

cv <- data %>% dplyr::filter(SET == "CV")

combine_cv <- combine_vars(data = cv) %>%
  dplyr::mutate(SET = "CV")

clin_var_levels <- c(
  "Overall",
  glue::glue("Age\nBinary"),
  "Sex",
  "Race",
  "Ethnicity",
  "Comorbidities",
  glue::glue("Prior\nCancer"),
  glue::glue("Smoking\nStatus"),
  glue::glue("Pack-\nYears")
)
category_levels <- c(
  "Overall",
  "Age >= 65 years",
  "Age < 65 years",
  "Female",
  "Male",
  "Other",
  "Black or African American",
  "White",
  " Other",
  "Hispanic or Latino",
  "Not Hispanic or Latino",
  "Overweight/Obesity, BMI > 25",
  "Type-II Diabetes",
  "COPD",
  "No Prior Cancer",
  "Prior Cancer",
  "Former",
  "Current",
  ">= 30 pack-years",
  "< 30 pack-years"
)

combine_cv_v5 <- combine_cv %>%
  dplyr::filter(
    !category %in% c("Unknown or not reported")
  ) %>%
  dplyr::mutate(
    clinical_var = factor(clinical_var, levels = clin_var_levels),
    category = factor(category, levels = category_levels),
  ) %>%
  dplyr::arrange(clinical_var, category)

ref_lines <- tibble::tibble(label = c("Specificity", "Sensitivity"), ref = c(0.5, 0.5))

ggplot(data = combine_cv_v5, aes(x = category, y = accuracy, ymin = ci_lower, ymax = ci_upper)) +
  geom_pointrange(
    aes(group = order, shape = label, color = label),
    size = point_size,
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(data = ref_lines, aes(yintercept = ref, color = label), linetype = "solid") +
  geom_hline(yintercept = c(0, 0.2, 0.4, 0.6, 0.8, 1), linetype = "dashed", linewidth = 0.2) +
  facet_grid(clinical_var ~ label, scales = "free", space = "free") +
  scale_x_discrete(
    breaks = category_levels,
    labels = stringr::str_wrap(paste0(category_levels, " (n=", combine_cv_v5$total[seq(1, nrow(combine_cv_v5), by = 2)], ")"), width = 20)
  ) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = font_size_x),
    strip.text.y = element_text(size = font_size_y),
    strip.background = element_blank(),
    axis.text.x = element_text(size = font_size_x),
    axis.text.y = element_text(size = font_size_y),
    panel.spacing = unit(1, "mm"),
    legend.position = "none"
  ) +
  scale_shape_manual(values = c(16, 16, 1)) +
  scale_color_manual(values = cbCol[c(6, 1)]) +
  labs(x = "", y = "")
