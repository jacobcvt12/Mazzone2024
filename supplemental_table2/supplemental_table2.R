library(dplyr)
library(forcats)
library(gt)
library(stringr)
library(tibble)

# Data is simulated for the purpose of code demonstration. It doesn't represent the real data.
data = readRDS("./data/cds_simulated_data.rds")

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

df <- data %>% 
  dplyr::filter(CNCLBL == "Lung Cancer") %>% 
  dplyr::mutate(
    SET = factor(SET, levels = c("TRAINING", "CV")),
    LCNCAT = factor(
      HISTOLOGY, 
      levels = c("Small Cell Carcinoma", "Non-Small Cell Carcinoma", "Unknown or not reported")
      ),
    LCNHIST_NEW = factor(
      LCNHIST_NEW,
      levels = c(
        "Adenocarcinoma", "Adenosquamous Carcinoma", "Large Cell Carcinoma",
        "Other", "Small Cell Carcinoma", "Squamous Cell Carcinoma", "Unknown or not reported"
        )
    ),
    LCNSTAGE = forcats::fct_drop(LCNSTAGE_NEW),
    LCNSTAGE2 = dplyr::if_else(
      is.na(`STAGE, IA OR IB-IV`), "Unknown or not reported", `STAGE, IA OR IB-IV`
    ),
    LCNSTAGE2 = factor(LCNSTAGE2, levels = c("IA", "IB-IV", "Unknown or not reported")),
    LCTDESC = factor(
      `TUMOR STAGE`, 
      levels = c("T1", "T2", "T3", "T4", "T0, TX, Unk"),
      labels = c("T1", "T2", "T3", "T4", "TX, T0, Unknown, or Not reported")
        ),
    LCNDESC = factor(
      `NODE STAGE`, 
      levels = c("N0", "N1", "N2", "N3", "NX, Unk"),
      labels = c("N0", "N1", "N2", "N3", "NX, Unknown, or Not reported")
      ),
    LCMDESC = factor(
      `METASTASIS STAGE`, 
      levels = c("M0", "M1", "MX, Unk"),
      labels = c("M0", "M1", "MX, Unknown, or Not reported"),
      )
  )

sub_tab2_lccat <- generate_cross_tab(df, "LCNCAT", "SET") %>% 
  dplyr::mutate(
    characteristic = dplyr::case_when(
      LCNCAT == "Small Cell Carcinoma" ~ "Histology",
      TRUE ~ " ")
  ) %>% 
  dplyr::select(characteristic,
                categories = LCNCAT, 
                stringr::str_subset(names(.), ".*\\(n.*\\)"))

sub_tab2_hist <- generate_cross_tab(df, "LCNHIST_NEW", "SET") %>% 
  dplyr::mutate(
    characteristic = dplyr::case_when(
      LCNHIST_NEW == "Adenocarcinoma" ~ "Histology Subtype",
      TRUE ~ " ")
  ) %>% 
  dplyr::select(characteristic,
                categories = LCNHIST_NEW, 
                stringr::str_subset(names(.), ".*\\(n.*\\)"))

sub_tab2_stg <- generate_cross_tab(df, "LCNSTAGE", "SET") %>% 
  dplyr::mutate(
    characteristic = dplyr::case_when(
      LCNSTAGE == "I" ~ "Stage",
      TRUE ~ " ")
  ) %>% 
  dplyr::select(characteristic,
                categories = LCNSTAGE, 
                stringr::str_subset(names(.), ".*\\(n.*\\)"))

sub_tab2_stg2 <- generate_cross_tab(df, "LCNSTAGE2", "SET") %>% 
  dplyr::mutate(
    characteristic = dplyr::case_when(
      LCNSTAGE2 == "IA" ~ "Stage, IA or IB-IV",
      TRUE ~ " ")
  ) %>% 
  dplyr::select(characteristic,
                categories = LCNSTAGE2, 
                stringr::str_subset(names(.), ".*\\(n.*\\)"))

sub_tab2_t <- generate_cross_tab(df, "LCTDESC", "SET") %>% 
  dplyr::mutate(
    characteristic = dplyr::case_when(
      LCTDESC == "T1" ~ "Tumor Stage",
      TRUE ~ " ")
  ) %>% 
  dplyr::select(characteristic,
                categories = LCTDESC, 
                stringr::str_subset(names(.), ".*\\(n.*\\)"))

sub_tab2_n <- generate_cross_tab(df, "LCNDESC", "SET") %>% 
  dplyr::mutate(
    characteristic = dplyr::case_when(
      LCNDESC == "N0" ~ "Node Stage",
      TRUE ~ " ")
  ) %>% 
  dplyr::select(characteristic,
                categories = LCNDESC, 
                stringr::str_subset(names(.), ".*\\(n.*\\)"))

sub_tab2_m <- generate_cross_tab(df, "LCMDESC", "SET") %>% 
  dplyr::mutate(
    characteristic = dplyr::case_when(
      LCMDESC == "M0" ~ "Metastasis Stage",
      TRUE ~ " ")
  ) %>% 
  dplyr::select(characteristic,
                categories = LCMDESC, 
                stringr::str_subset(names(.), ".*\\(n.*\\)"))

sub_tab2_df <- list(sub_tab2_stg, sub_tab2_stg2, sub_tab2_t, sub_tab2_n, sub_tab2_m,
                sub_tab2_lccat, sub_tab2_hist) %>% 
  purrr::reduce(dplyr::bind_rows) %>% 
  dplyr::rename(" " = characteristic, 
                "  " = categories)

sup_tab2 <- gt::gt(sub_tab2_df)
sup_tab2
