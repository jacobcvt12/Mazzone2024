library(dplyr)
library(forcats)
library(ggplot2)

# All fields included in this dataset are simulated for the purpose of code demonstration 
# and do not represent real results or actual data collected participants.
data = readRDS("./data/cds_simulated_data.rds")
cv <- data %>% dplyr::filter(SET == "CV")

font_size <- 5

ref <- 0.22
cv_stage <- cv %>%
  dplyr::filter(LCNSTAGE_NEW != "Unknown or not reported") %>%
  dplyr::group_by(LCNSTAGE_NEW) %>%
  dplyr::mutate(LCNSTAGE_NEW_CT = paste0(LCNSTAGE_NEW, "\n", "(n=", dplyr::n(), ")")) %>%
  dplyr::arrange(LCNSTAGE_NEW) %>%
  dplyr::mutate(LCNSTAGE_NEW_CT = forcats::fct_inorder(LCNSTAGE_NEW_CT))

ggplot(cv_stage, aes(x = LCNSTAGE_NEW_CT, y = DELFI_SCORE)) +
  geom_boxplot(alpha = 0) + 
  geom_jitter(position = position_jitter(0.3), alpha = 0.2) +
  geom_hline(yintercept = ref, linetype = "dashed") +
  scale_x_discrete(name = "") +
  labs(y = "Delfi Score") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(
    axis.title = element_text(size = font_size),
    axis.text.x = element_text(size = font_size),
    axis.text.y = element_text(size = font_size)
  )


cv_hist <- cv %>%
  dplyr::mutate(
    LCNHIST_NEW = forcats::fct_collapse(
      LCNHIST_NEW,
      "Non-cancer" = "Non-Cancer",
      "Adeno" = "Adenocarcinoma",
      "Squamous" = "Squamous Cell Carcinoma",
      "Small cell" = "Small Cell Carcinoma",
      "Large cell" = "Large Cell Carcinoma",
      "Adenosquamous" = "Adenosquamous Carcinoma",
      "Other" = c("Other", "Unknown or not reported")
    )
  ) %>%
  dplyr::group_by(LCNHIST_NEW) %>%
  dplyr::mutate(LCNHIST_NEW_CT = paste0(LCNHIST_NEW, "\n", "(n=", dplyr::n(), ")")) %>%
  dplyr::arrange(LCNHIST_NEW) %>%
  dplyr::mutate(LCNHIST_NEW_CT = forcats::fct_inorder(LCNHIST_NEW_CT))


ggplot(cv_hist, aes(x = LCNHIST_NEW_CT, y = DELFI_SCORE)) +
  geom_boxplot(alpha = 0) + 
  geom_jitter(position = position_jitter(0.3), alpha = 0.2) +
  geom_hline(yintercept = ref, linetype = "dashed") +
  labs(x = "", y = "Delfi Score") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(
    axis.title = element_text(size = font_size),
    axis.text.x = element_text(size = font_size),
    axis.text.y = element_text(size = font_size)
  )
