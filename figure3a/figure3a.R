library(readxl)
library(tidyverse)
library(recipes)
library(circlize)
library(ComplexHeatmap)
library(stringr)

dat.scores <- 
  read_csv("./data/prediction_lucas.csv") %>% select(id, `Delfi Score` = score.seq)

dat.meta <- 
  read_excel("./data/LUCAS_metadata.xlsx", sheet = 1) %>% 
  filter(id %in% dat.scores$id) %>%
  rename(status = type, stage = Stage)

dat.labels <- 
  dat.meta %>% 
  select(id, status, stage) 

dat.features.raw <- 
  read_csv("./data/training-set.csv") 

dat.full <- 
  dat.features.raw %>% 
  inner_join(dat.labels) %>% 
  select(id, status, stage, 
         starts_with("ratio_"), starts_with("zscore_"), 
         Age = clinical_age, `Smoking Status` = clinical_smokingstatus, `Pack Years` = clinical_packyears) %>%
  mutate(`Pack Years` = case_when(`Pack Years` <= 30 ~ "30 or fewer",
                                  `Pack Years` >= 31 & `Pack Years` <= 60 ~ "31 - 60",
                                  `Pack Years` >= 61 & `Pack Years` <= 90 ~ "61 - 90",
                                  `Pack Years` >= 91 & `Pack Years` <= 120 ~ "91 - 120",
                                  `Pack Years` >= 121 ~ "121 or more")) %>% 
  mutate(`Pack Years` = factor(`Pack Years`, levels = rev(c("30 or fewer", "31 - 60", "61 - 90", "91 - 120", "121 or more")))) %>% 
  mutate(`Smoking Status` = str_to_title(`Smoking Status`)) %>% 
  mutate(`Smoking Status` = factor(`Smoking Status`, levels = c("Never", "Former", "Current")))

dat.features.all <- 
  dat.full %>% 
  select(id, status, starts_with("ratio_"), starts_with("zscore_"))

dat.features.noncancer <- 
  dat.full %>% 
  filter(status == "healthy") %>% 
  select(id, status, starts_with("ratio_"), starts_with("zscore_"))

dat.features.noncancer <- 
  dat.full %>% 
  filter(status == "healthy") %>% 
  select(id, status, starts_with("ratio_"), starts_with("zscore_"))

dat.clinical <- 
  dat.full %>% 
  select(id, Age, `Smoking Status`, `Pack Years`)

recipe.noncancer <- 
  recipe(status ~ ., data = dat.features.noncancer) %>%
  update_role(id, new_role="ID") %>%
  step_normalize(all_predictors()) %>%
  prep()

baked.full <-
  recipe.noncancer %>%
  bake(new_data = dat.features.all) %>%
  mutate(id = dat.features.all$id) %>%
  mutate(status = dat.features.all$status) %>% 
  inner_join(dat.scores) %>%
  arrange(desc(`Delfi Score`)) %>%
  inner_join(dat.meta) %>% 
  mutate(Status = if_else(status == "cancer", "Lung Cancer", "Non-cancer")) %>%
  mutate(Status = factor(Status, levels = c("Non-cancer", "Lung Cancer"))) %>% 
  mutate(Stage = case_when(stage %in% c("IA", "IB") ~ "I",
                           stage %in% c("IIA", "IIB") ~ "II",
                           stage %in% c("IIIA", "IIIB", "IIIC") ~ "III",
                           stage %in% c("IV") ~ "IV",
                           TRUE ~ "Non-cancer")) %>% 
  inner_join(dat.clinical)

baked.heatmap <- 
  baked.full %>% 
  select(starts_with("ratio_"), starts_with("zscore_")) %>% 
  as.matrix()

annot.rows <- 
  baked.full[, c('Delfi Score', 'Status', 'Stage', 'Age', 'Smoking Status', 'Pack Years')]

col.fun.row = list(Age = colorRamp2(c(23, 94), c('#f7fcf0', '#4eb3d3')) ,
                   `Smoking Status` = structure(c('#ffffff', '#737373', 'grey20'), names = c("Never", "Former", "Current")),
                   `Pack Years` = structure(colorRampPalette(c("#f6eff7", "#3690c0"))(5), names = c("30 or fewer", "31 - 60", "61 - 90", "91 - 120", "121 or more")),
                   Status = structure(c('#f7f7f7', '#238b45'), names = c("Non-cancer", "Lung Cancer")),
                   Stage = structure(c('#eeeeee', '#fee5d9', '#fcae91', '#fb6a4a', '#cb181d'), names = c('Non-cancer', 'I', 'II','III','IV')),
                   `Delfi Score` = colorRamp2(c(0, 1), c('#ffffff', '#111111')))

annot.cols <-
  data.frame(Family = sapply(colnames(baked.heatmap),
                             function(x) {
                               
                               case_when(str_detect(x, "ratio_") ~ "Fragmentation profile",
                                         str_detect(x, "zscore_") ~ "Chromosomal changes")
                               
                             })) %>%
  mutate(Family = factor(Family, levels = c("Fragmentation profile", "Chromosomal changes")))

names(annot.cols) <- "Feature Classes"
col.fun.column = list(`Feature Classes` = structure(c("#8da0cb", "#66c2a5"), names = c("Fragmentation profile", "Chromosomal changes")))

col.fun.body = colorRamp2(c(-2, -1, 0, 1, 2), rev(c('#d7191c', '#fdae61', '#ffffbf', '#abd9e9', '#2c7bb6')))

ht.raw <- 
  Heatmap(matrix = baked.heatmap, 
          col = col.fun.body,
          name = "Values", 
          row_title = NULL,
          column_title = NULL,
          show_row_dend =  FALSE, 
          show_column_dend = FALSE,
          show_row_names = FALSE, 
          show_column_names = FALSE,
          left_annotation = rowAnnotation(df = annot.rows, 
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

ht.out <- 
  columnAnnotation(df = annot.cols, 
                   col = col.fun.column, 
                   annotation_name_side = 'right',
                   annotation_name_gp = gpar(fontsize = 7),
                   annotation_legend_param = list(labels_gp = gpar(fontsize = 7),
                                                  title_gp = gpar(fontsize = 7,
                                                                  fontface = "bold"))) %v% 
  ht.raw

draw(ht.out, merge_legends = TRUE, padding = unit(c(1, 1, 27, 1), "mm"))
pdf("figure3a.pdf", height=6.69291, width=7.08661)
