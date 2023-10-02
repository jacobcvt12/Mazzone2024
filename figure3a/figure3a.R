library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(recipes)

dat.ids <- read_csv("training_ids.csv") %>% select(id, USUBJID)
dat.ct <- read_csv("dat_ct.csv")
dat.ct.noncancer <- dat.ct %>% filter(status == 0)
dat.cv <- read_csv("dat_cv.csv")

dat.meta <- 
  read_csv("./combined_ct_cv_metadata.csv") %>% 
  inner_join(dat.ids) %>% 
  select(id, Stage = LCNSTAGE, Age = AGE, Sex = SEX, `Pack Years` = SUPACKYR, `Smoking Status` = SUCAT, Histology=LCNHIST) %>%
  mutate(`Pack Years` = case_when(`Pack Years` <= 30 ~ "30 or fewer",
                                  `Pack Years` >= 31 & `Pack Years` <= 60 ~ "31 - 60",
                                  `Pack Years` >= 61 & `Pack Years` <= 90 ~ "61 - 90",
                                  `Pack Years` >= 91 & `Pack Years` <= 120 ~ "91 - 120",
                                  `Pack Years` >= 121 ~ "121 or more"),
         Histology = case_when(Histology == "Adenocarcinoma" ~ "Adenocarcinoma",
                               Histology == "Squamous cell carcinoma" ~ "Squamous cell carcinoma",
                               Histology == "Not Confirmed Lung Cancer" ~ "Non-cancer",
                               TRUE ~ "Other")) %>%
  mutate(`Pack Years` = factor(`Pack Years`, levels = rev(c("30 or fewer", "31 - 60", "61 - 90", "91 - 120", "121 or more"))))

scores <- 
  read_csv("final_scores.csv") %>% 
  select(id, `Delfi Score` = score) %>%
  inner_join(dat.meta)

recipe.ct.noncancer <- 
  recipe(status ~ ., data = dat.ct.noncancer) %>%
  update_role(id, new_role="ID") %>%
  step_normalize(all_predictors()) %>%
  prep()

baked.ct <-
  recipe.ct.noncancer %>%
  bake(new_data = dat.ct) %>%
  mutate(set = "CT") %>%
  mutate(id = dat.ct$id) %>%
  inner_join(scores) %>%
  arrange(desc(`Delfi Score`)) %>%
  mutate(Status = if_else(status == 1, "Lung Cancer", "Non-cancer")) %>%
  mutate(Status = factor(Status, levels = c("Non-cancer", "Lung Cancer"))) %>%
  mutate(Stage = if_else(Stage == "Not Confirmed Lung Cancer", "Non-cancer", Stage)) 

baked.ct.heatmap <- 
  baked.ct %>% 
  select(contains(c("centeredslratio", "zscore", "mxab", "mtDNA", "mixture"))) %>% 
  as.matrix()

annot.rows.ct <- 
  baked.ct[, c('Delfi Score', 'Status', 'Stage', 'Histology', 'Age', 'Sex', 'Smoking Status', 'Pack Years')]

col.fun.row = list(Age = colorRamp2(c(50, 90), c('#f7fcf0', '#4eb3d3')) ,
                   Sex = structure(c('#807dba', '#dadaeb'), names = c('Male', 'Female')),
                   Histology = structure(c('#faebd7', '#ee82ee', '#483d8b', '#add8e6'), names=c('Non-cancer', 'Adenocarcinoma', 'Squamous cell carcinoma', 'Other')),
                   `Smoking Status` = structure(c('#ffffff', '#737373'), names = c("Former", "Current")),
                   `Pack Years` = structure(colorRampPalette(c("#f6eff7", "#3690c0"))(5), names = c("30 or fewer", "31 - 60", "61 - 90", "91 - 120", "121 or more")),
                   Status = structure(c('#f7f7f7', '#238b45'), names = c("Non-cancer", "Lung Cancer")),
                   Stage = structure(c('#eeeeee', '#fee5d9', '#fcae91', '#fb6a4a', '#cb181d'), names = c('Non-cancer', 'I', 'II','III','IV')),
                   `Delfi Score` = colorRamp2(c(0, 1), c('#ffffff', '#111111')))
col.fun.column = list(Family = structure(c('#66c2a5', '#8da0cb', "#D55E00", "#CC79A7"), names = c('Chromosomal changes', 'Fragmentation profile', "mtDNA", "Length distribution")))
col.fun.body = colorRamp2(c(-2, -1, 0, 1, 2), rev(c('#d7191c', '#fdae61', '#ffffbf', '#abd9e9', '#2c7bb6')))

ht.raw <- 
  Heatmap(matrix = baked.ct.heatmap, 
          col = col.fun.body,
          name = "Values", 
          row_title = NULL,
          column_title = NULL,
          show_row_dend =  FALSE, 
          show_column_dend = FALSE,
          show_row_names = FALSE, 
          show_column_names = FALSE,
          left_annotation = rowAnnotation(df = annot.rows.ct, 
                                          col = col.fun.row, 
                                          annotation_name_side= 'top', 
                                          annotation_name_offset = unit(3, "mm"),
                                          annotation_name_gp = gpar(fontsize = 7),
                                          annotation_legend_param = list(labels_gp = gpar(fontsize = 7),
                                                                         title_gp = gpar(fontsize = 7,
                                                                                         fontface = "bold"))),
          cluster_row_slices = FALSE, 
          cluster_rows = FALSE,
          cluster_column_slices = FALSE,
          cluster_columns = FALSE,
          gap = unit(2, "mm"),
          column_names_max_height = unit(5, "mm"),
          heatmap_legend_param = list(labels_gp = gpar(fontsize = 7),
                                      title_gp = gpar(fontsize = 7,
                                                      fontface = "bold")))

ht <- 
  columnAnnotation(df = annot.cols, 
                   col = col.fun.column, 
                   annotation_name_side = 'right',
                   annotation_name_gp = gpar(fontsize = 7),
                   annotation_legend_param = list(labels_gp = gpar(fontsize = 7),
                                                  title_gp = gpar(fontsize = 7,
                                                                  fontface = "bold"))) %v% 
  ht.raw

draw(ht, merge_legends = TRUE, padding = unit(c(1, 1, 27, 1), "mm"))
pdf("figure3a.pdf", height=6.69291, width=7.08661)

