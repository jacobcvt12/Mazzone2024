library(dplyr)
library(gt)
library(purrr)
library(rlang)
library(stringr)
library(tibble)
library(tidyr)

# All fields included in this dataset are simulated for the purpose of code demonstration 
# and do not represent real results or actual data collected from participants.
data = readRDS("./data/cds_simulated_data.rds")

data <- data %>%
  dplyr::mutate(SET = factor(SET, levels = c("TRAINING", "CV")))

format_tab1 <- function(df, cat_var, pattern, val, lab) {
  df %>%
    dplyr::select(
      categories = all_of(cat_var),
      stringr::str_subset(names(.), pattern),
      Total
    ) %>%
    dplyr::mutate(
      characteristic = dplyr::if_else(categories == val, lab, ""),
      .before = 1
    )
}

calculate_quantiles <- function(df, num_var, cat_var, digit = 2) {
  if (cat_var == "") {
    quants <- df %>% 
      dplyr::summarise(
        med = round(quantile(!!as.symbol(num_var), probs = 0.5, na.rm = TRUE), digit),
        q1 = round(quantile(!!as.symbol(num_var), probs = 0.25, na.rm = TRUE), digit),
        q3 = round(quantile(!!as.symbol(num_var), probs = 0.75, na.rm = TRUE), digit),
        total = n()
      ) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(
        lab = "Total",
        value = paste0(med, " (", q1, ", ", q3, ")")
      )
  } else {
    quants <- df %>% 
      dplyr::group_by(!!as.symbol(cat_var)) %>% 
      dplyr::summarise(
        med = round(quantile(!!as.symbol(num_var), probs = 0.5, na.rm = TRUE), digit),
        q1 = round(quantile(!!as.symbol(num_var), probs = 0.25, na.rm = TRUE), digit),
        q3 = round(quantile(!!as.symbol(num_var), probs = 0.75, na.rm = TRUE), digit),
        total = n()
      ) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(
        lab = paste0(!!as.symbol(cat_var), " (n=", total, ")"),
        value = paste0(med, " (", q1, ", ", q3, ")")
      ) 
  }
  
  quants <- quants %>%   
    dplyr::select(lab, value) %>% 
    tidyr::pivot_wider(names_from = lab, values_from = value) 
  
  return(quants)
}

calculate_row_total <- function(df, category1 = "Lung Cancer", category2 = "Non-Cancer") {
  df %>% 
    dplyr::mutate(total = sum(!!rlang::sym(category1)) + sum(!!rlang::sym(category2))) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(
      row_total = !!rlang::sym(category1) + !!rlang::sym(category2),
      row_total_perc = round(row_total / total * 100, 0),
      Total = paste0(row_total, " (", row_total_perc, "%)")
    ) %>% 
    dplyr::ungroup()
}

generate_cross_tab <- function(df, cat_var1, cat_var2) {
  cat_var1_lvls <- tibble::tibble(!!as.symbol(cat_var1) := levels(df[[cat_var1]]))
  
  tab <- cat_var1_lvls %>% 
    dplyr::left_join(
      df %>% dplyr::count(!!as.symbol(cat_var1), !!as.symbol(cat_var2)),
      by = cat_var1) %>%
    dplyr::filter(!(!!as.symbol(cat_var1) == "Unknown or not reported" & is.na(n))) %>% 
    tidyr::pivot_wider(names_from = !!as.symbol(cat_var2), values_from = n) %>% 
    dplyr::select(-matches("NA"))
  
  tab[is.na(tab)] = 0
  
  for (i in 2:length(tab)) {
    new_col <- paste0(names(tab)[i], " (n=", sum(tab[i]), ")")
    
    tab <- tab %>% 
      dplyr::mutate(
        !!paste0("cat", i-1, "_total") := sum(!!as.symbol(names(.)[i]))
      ) %>%
      dplyr::rowwise() %>% 
      dplyr::mutate(
        !!paste0("cat", i-1, "_perc") := round(
          !!as.symbol(names(.)[i]) / !!as.symbol(paste0("cat", i-1, "_total")), 2
        ) * 100
        ,
        !!new_col := paste0(
          !!as.symbol(names(.)[i]), " (", !!as.symbol(paste0("cat", i-1, "_perc")), "%)"
        )
      ) %>% 
      dplyr::ungroup()
  }
  
  tab <- tab %>% 
    dplyr::select(-c(ends_with("total")))
  
  return(tab)
}

all_tab1 <- data %>%
  dplyr::mutate(
    CNCLBL = dplyr::if_else(CNCLBL == "Non-cancer", "Non-Cancer", CNCLBL)
  )

tab1_age_num <- calculate_quantiles(df = all_tab1, num_var = "AGE", cat_var = "CNCLBL") %>%
  dplyr::bind_cols(calculate_quantiles(df = all_tab1, num_var = "AGE", cat_var = "")) %>%
  dplyr::mutate(characteristic = "Age (years), median (Q1, Q3)", categories = "", .before = 1)

tab1_agecat <- generate_cross_tab(all_tab1, "AGECAT", "CNCLBL") %>%
  calculate_row_total() %>%
  format_tab1(
    cat_var = "AGECAT",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "Age < 65 years",
    lab = "Age, N (%)"
  )

tab1_sex <- generate_cross_tab(all_tab1, "SEX", "CNCLBL") %>%
  calculate_row_total() %>%
  format_tab1(
    cat_var = "SEX",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "Male",
    lab = "Sex, N (%)"
  )

tab1_race <- generate_cross_tab(all_tab1, "RACE_NEW", "CNCLBL") %>%
  calculate_row_total() %>%
  format_tab1(
    cat_var = "RACE_NEW",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "White",
    lab = "Race, N (%)"
  )

tab1_ethnic <- generate_cross_tab(all_tab1, "ETHNIC", "CNCLBL") %>%
  calculate_row_total() %>%
  format_tab1(
    cat_var = "ETHNIC",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "Not Hispanic or Latino",
    lab = "Ethnicity, N (%)"
  )

tab1_edu <- generate_cross_tab(all_tab1, "EDU", "CNCLBL") %>%
  calculate_row_total() %>%
  format_tab1(
    cat_var = "EDU",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "Less than High School Graduate",
    lab = "Education, N (%)"
  )

tab1_region <- generate_cross_tab(all_tab1, "REGION", "CNCLBL") %>%
  calculate_row_total(category1 = "Lung Cancer", category2 = "Non-Cancer") %>%
  format_tab1(
    cat_var = "REGION",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "Northeast",
    lab = "Geographic Region"
  )

tab1_bmi_num <- calculate_quantiles(df = all_tab1, num_var = "BMI", cat_var = "CNCLBL") %>%
  dplyr::bind_cols(calculate_quantiles(df = all_tab1, num_var = "BMI", cat_var = "")) %>%
  dplyr::mutate(characteristic = "BMI, median (Q1, Q3)", categories = "", .before = 1)

tab1_new_uspstf <- generate_cross_tab(all_tab1, "USPSTFNC", "CNCLBL") %>%
  calculate_row_total() %>%
  format_tab1(
    cat_var = "USPSTFNC",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "USPSTF (2021)",
    lab = "Lung Cancer USPSTF 2021 Criteria, N (%)"
  )

tab1_pcyr_num <- calculate_quantiles(df = all_tab1, num_var = "SUPACKYR", cat_var = "CNCLBL") %>%
  dplyr::bind_cols(calculate_quantiles(df = all_tab1, num_var = "SUPACKYR", cat_var = "")) %>%
  dplyr::mutate(characteristic = "Pack-Years, median (Q1, Q3)", categories = "", .before = 1)

tab1_smkint_num <- calculate_quantiles(df = all_tab1, num_var = "CIG_PER_DAY", cat_var = "CNCLBL") %>%
  dplyr::bind_cols(calculate_quantiles(df = all_tab1, num_var = "CIG_PER_DAY", cat_var = "")) %>%
  dplyr::mutate(
    characteristic = "Smoking Cigarettes per Day, median (Q1, Q3)",
    categories = "", .before = 1
  )

tab1_quit_num <- calculate_quantiles(df = all_tab1, num_var = "SUQUIT", cat_var = "CNCLBL") %>%
  dplyr::bind_cols(calculate_quantiles(df = all_tab1, num_var = "SUQUIT", cat_var = "")) %>%
  dplyr::mutate(
    characteristic = "Years Since Smoking Cessation, median (Q1, Q3)",
    categories = "", .before = 1
  )

tab1_smkstat <- generate_cross_tab(all_tab1, "SUCAT", "CNCLBL") %>%
  calculate_row_total() %>%
  format_tab1(
    cat_var = "SUCAT",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "Current",
    lab = "Smoking Status, N (%)"
  )

tab1_copd <- generate_cross_tab(all_tab1, "MH_COPD", "CNCLBL") %>%
  calculate_row_total() %>%
  format_tab1(
    cat_var = "MH_COPD",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "COPD",
    lab = "COPD Status, N (%)"
  )

tab1_diab <- generate_cross_tab(all_tab1, "MH_DIAB", "CNCLBL") %>%
  calculate_row_total() %>%
  format_tab1(
    cat_var = "MH_DIAB",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "Type-II Diabetes",
    lab = "Diabetes Status, N (%)"
  )

tab1_prior_nlc <- generate_cross_tab(all_tab1, "PROTHCYN", "CNCLBL") %>%
  calculate_row_total() %>%
  format_tab1(
    cat_var = "PROTHCYN",
    pattern = "(-)?Cancer \\(.*\\)",
    val = "Yes",
    lab = "History of Non-Lung Cancer, N (%)"
  )

all_tab1_df <- list(
  tab1_age_num, tab1_agecat, tab1_sex,
  tab1_race, tab1_ethnic, tab1_edu, 
  tab1_region, tab1_bmi_num, tab1_new_uspstf,
  tab1_pcyr_num, tab1_smkint_num, tab1_quit_num, 
  tab1_smkstat, tab1_copd, tab1_diab, 
  tab1_prior_nlc
) %>%
  purrr::reduce(dplyr::bind_rows) %>%
  dplyr::mutate(!!paste0("Total (n=", nrow(all_tab1), ")") := Total) %>%
  dplyr::select(-Total) %>%
  dplyr::rename(" " = characteristic, "  " = categories)

sup_tab1 <- gt::gt(all_tab1_df)

sup_tab1
