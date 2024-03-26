library(dplyr)
library(gt)
library(purrr)
library(rlang)
library(stringr)
library(tibble)
library(tidyr)

# Data is simulated for the purpose of code demonstration. It doesn't represent the real data.
data = readRDS("./data/cds_simulated_data.rds")

data <- data %>%
  dplyr::mutate(SET = factor(SET, levels = c("TRAINING", "CV"))) %>%
  dplyr::mutate(
    nod_6mm = dplyr::case_when(
      ENRLLNLA >= 6 ~ ">= 6",
      ENRLLNLA < 6 ~ "< 6",
      TRUE ~ "Unknown"
    ),
    rad_nod = dplyr::case_when(
      stringr::str_detect(RADSORRESBC, "1") ~ "Lung-RADS 1 (No nodules or nodules with specific calcifications)",
      stringr::str_detect(RADSORRESBC, "2") ~ "Lung-RADS 2 or nodule(s) <6 mm if no Lung-RADS",
      stringr::str_detect(RADSORRESBC, "3|4") ~ "Lung-RADS 3-4 or nodule(s) >=6mm if no Lung-RADS",
      nod_6mm == "< 6" ~ "Lung-RADS 2 or nodule(s) <6 mm if no Lung-RADS",
      nod_6mm == ">= 6" ~ "Lung-RADS 3-4 or nodule(s) >=6mm if no Lung-RADS",
      TRUE ~ "Lung-RADS and nodule size unknown or not reported"
    ),
    rad_nod = factor(rad_nod, levels = c(
      "Lung-RADS 1 (No nodules or nodules with specific calcifications)",
      "Lung-RADS 2 or nodule(s) <6 mm if no Lung-RADS",
      "Lung-RADS 3-4 or nodule(s) >=6mm if no Lung-RADS",
      "Lung-RADS and nodule size unknown or not reported"
      ))
    )

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

calculate_classifier_accuracy <- function(df, pred = c("pos", "neg"), type = "sens") {
  if (type == "sens") {
    prediction <- stringr::str_subset(pred, "(P|p)ositive")
    new_col <- "sens_CIs"
  } else if (type == "spec"){
    prediction <- stringr::str_subset(pred, "(N|n)egative")
    new_col <- "spec_CIs"
  } else {
    stop("Type argument can only be 'sens' or 'spec'.")
  }
  
  df <- df %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(
      total = sum(!!as.symbol(pred[1]), !!as.symbol(pred[2])),
      pred_accuracy = round(!!as.symbol(prediction) / total, 2),
      lower = binom::binom.wilson(x = !!as.symbol(prediction), n = total, conf.level = 0.95)["lower"],
      upper = binom::binom.wilson(x = !!as.symbol(prediction), n = total, conf.level = 0.95)["upper"],
      !!(new_col) := paste0(
        sprintf("%.2f", pred_accuracy), " (", sprintf("%.2f", lower), ", ", sprintf("%.2f", upper), ")"
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-c(total, pred_accuracy, lower, upper))
  
  return(df)
}


tab3_df <- data %>% 
  dplyr::filter(SET == "CV", CNCLBL == "Non-cancer") %>% 
  generate_cross_tab("rad_nod", "DELFI_RESULT") %>% 
  calculate_classifier_accuracy(pred = c("DLCST Positive", "DLCST Negative"), type = "spec") %>% 
  dplyr::select("Enrollment CT Scan Result" = rad_nod, 
                stringr::str_subset(names(.), "Positive \\(n.*\\)"),
                stringr::str_subset(names(.), "Negative \\(n.*\\)"),
                "Specificity (95% CI)" = spec_CIs) 
names(tab3_df) <- stringr::str_replace(names(tab3_df), "DLCST", "Test")

tab3 <- gt::gt(tab3_df)
tab3
