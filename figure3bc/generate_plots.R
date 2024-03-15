
library(RCurl)
library(recipes)
library(reshape2)
library(tidyverse)
library(readxl)
library(caret)

####################
# Load data and model
####################

#extract features
#read in model object
trained.model <- readRDS(url("https://github.com/cancer-genomics/reproduce_lucas_wflow/raw/master/code/models_c1/model_seq_glm.rds", method="libcurl"))
features <- trained.model$trainingData %>%
  mutate(id=factor(id),
         type=factor(type))

#set recipe
recipe_seq <- recipe(type ~ ., data=features) %>%
  step_rm(starts_with("clinical_"), multinucratio) %>%
  update_role(id, new_role = "ID") %>%
  step_pca(starts_with("ratio_"), prefix = "ratio_pc_",  threshold=0.90)  %>%
  step_corr(all_predictors(), threshold=0.95) %>%
  step_nzv(all_predictors())

#load prepped recipe
prepped.recipe<- prep(recipe_seq)
baked.features <- bake(prepped.recipe, new_data = features)

glmnetGrid <- expand.grid(
  alpha = 1,
  lambda = 10^seq(-5, -1, length.out = 100))
#### Train models
set.seed(1234)
ctrl <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 10,
                     verboseIter = FALSE,
                     savePredictions="final",
                     classProbs=TRUE,
                     index=createMultiFolds(features$type, 5, 10),
                     summaryFunction = twoClassSummary)

model_seq <- caret::train(recipe_seq,
                          data = features,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl,
                          metric="ROC")


#pull out coefficients
trained.model$levels
# [1] "cancer"  "healthy"
# this indicates that the higher level is healthy and that's coded as 1 in logistic regression
# need to flip the coefficient signs
model.coefficients <- coef(model_seq$finalModel, s = model_seq$bestTune$lambda) * (-1)

feature_means <- baked.features  %>%
  select(-c(id, type)) %>%
  colMeans()
feature_sds <- baked.features %>%
  select(-c(id, type)) %>%
  as.data.frame() %>%
  summarise_all(sd)

coefficients <- data.frame(features = names(feature_sds),
                            sd = as.numeric(feature_sds))
coefficients <- merge(coefficients,
                       data.frame(features = rownames(model.coefficients),
                                  model.coefficients = as.numeric(model.coefficients)),
                       by = 'features', all.x = TRUE)

coefficients$scaled_coefficients <- coefficients$model.coefficients * coefficients$sd

# PC heatmap-----------
# extract the PC loadings
loadings <- prepped.recipe %>% tidy(number = 2)
loadings <- data.frame(loadings)
loadings <- dcast(loadings, terms ~ component, value.var = 'value')
loadings$bin.id <- as.numeric(gsub('ratio_', '', loadings$terms))
loadings <- loadings[with(loadings, order(bin.id)),]
# only keep PC1 to PC11 given that these terms make the cut for the threshold
n.pcs.retained <- (prepped.recipe %>% tidy(number = 2, type = "variance") %>% filter(terms == "percent variance") %>% mutate(total_variance = cumsum(value)) %>% filter(total_variance > 90) %>% pull(component))[1]

loadings <- loadings[,c('terms', 'bin.id', sapply(seq(n.pcs.retained), function(x) paste0('PC', x)))]


# now work on visualization of the loadings
pd <- loadings[,! colnames(loadings) %in% c('terms')]
pd <- melt(pd, id.vars = c('bin.id'))
pd$pc.id <- as.numeric(gsub('PC', '', pd$variable))

# annotate pd with chromosome arm and position
bins5mb <- read.csv(text = getURL("https://raw.githubusercontent.com/cancer-genomics/reproduce_lucas_wflow/master/data/lucas_5mbs_delfi473.csv"))
bins5mb <- bins5mb %>% select(chr, start, end, bin, arm) %>% unique() %>% 
  filter(chr != "chrX")
locs <- bins5mb
locs$pos <- apply(locs[,c('start', 'end')], 1, mean)
pd <- merge(pd, locs[,c('chr','pos','arm', 'bin')],
            by.x = 'bin.id', by.y = 'bin', all.x = TRUE)
pd$arm <- factor(pd$arm, levels = unique(locs$arm))


# Importance barplot------------
sc <- coefficients %>% select(features, scaled_coefficients)
sc$abs.value = abs(sc$scaled_coefficients)
sc$sign.value = factor(sign(sc$scaled_coefficients), levels = c(-1, 1))
sc$feature.type <- sapply(sc$feature, function(x) strsplit(as.character(x), split = '_')[[1]][1])
sc$feature.type <- factor(sc$feature.type, levels = c('zscore', 'ratio'))
sc <- sc[with(sc, order(-abs.value, feature.type)),]
sc$features <- factor(sc$features, levels = sc$features)
sc <- subset(sc, scaled_coefficients != 0)

# determine feature order for the heatmap figure
f <- levels(sc$features)
pd$pc.id <- factor(as.character(pd$pc.id),
                   levels = as.character(as.numeric(gsub('ratio_pc_', '', f[grepl('ratio', f)]))))
pd$variable <- factor(as.character(pd$variable),
                   levels = as.character(paste0("PC", as.numeric(gsub('ratio_pc_', '', f[grepl('ratio', f)])))))

sc$features <- gsub('zscore_', 'Z ', gsub('ratio_pc_', 'PC ', sc$features))
sc$features <- gsub('PC 0', 'PC', sc$features)
sc$features <- factor(sc$features , levels = sc$features)

#---------------- now make plots -----------------#
# visualize feature by PC heatmap
heatmap <- ggplot(pd, aes(x = pc.id, y = bin.id)) +
  facet_grid(arm ~ variable, scales = 'free', switch = 'y')+
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = 'RdBu') +
  scale_y_continuous(trans = 'reverse') +
  theme_minimal() +
  labs(x = 'Principal components\nof fragmentation profiles', y = '', fill = 'Value') +
  theme(strip.text.x = element_text(size = 5, angle = 0, hjust = 0.5, vjust = 0.5),
        strip.text.y= element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid =  element_blank(),
        panel.spacing=unit(0.1, "lines"),
        legend.position = 'none') + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

#-------------------------------------------------------------------------------------------------------#
# top annotation
ta <- subset(sc, grepl('PC', features))
ta$pc.id <- gsub('PC', '', gsub(' ','',ta$features))
ta$sign.value = factor(ta$sign.value, levels = c(-1, 1))
ta <- ta[,c('pc.id', 'abs.value', 'sign.value')]
ta <- rbind(ta, data.frame(pc.id = setdiff(seq(n.pcs.retained), as.numeric(ta$pc.id)), abs.value =0, sign.value = NA))
ta$pc.id <- factor(ta$pc.id, levels = ta$pc.id)
ta$x = 1

cols <- c('#0571b0', '#ca0020', '#999999')
names(cols) <- c('-1', '1', NA)
ta$arm = 1
top.annot <- ggplot(ta, aes(x = pc.id, y = abs.value, colour = sign.value)) +
  facet_grid(arm ~ pc.id, scales = 'free', switch = 'y', space = 'free') +
  geom_point() +
  geom_segment(aes(x = pc.id, xend = pc.id, y = abs.value, yend = 0), linewidth = 1) +
  theme_minimal() +
  labs(x = '', y = 'Log odds ratio', fill = '') +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.text.y=element_blank(),
        strip.text.x = element_blank(),
        axis.line.x=element_line(color="gray"),
        axis.ticks.x=element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor.y  = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.spacing=unit(0.1, "lines"),
        legend.position = 'none') +
  scale_color_manual(values = cols, na.value = '#999999') + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#-------------------------------------------------------------------------------------------------------#
# side annotation --> arm importance
bl <- unique(pd[,c('chr','pos', 'bin.id', 'arm')])
bl$x = 1
arm.imp <- bl %>% group_by(arm) %>% 
  summarize(bin.id = mean(bin.id), x = mean(x))

arm.coef <- subset(sc, ! grepl('PC', features))[,c('features', 'abs.value', 'sign.value')]
arm.coef$arm <- gsub('Z ', '', arm.coef$features)
arm.coef$features <- NULL
arm.coef <- rbind(arm.coef,
                  data.frame(abs.value = 0, sign.value = NA, arm = setdiff(unique(locs$arm), arm.coef$arm)))

arm.imp <- merge(arm.imp, arm.coef[,c('arm', 'abs.value', 'sign.value')], by = 'arm', all.x = TRUE)

bl$pc.id <- 1
right.annot <- ggplot() +
  geom_point(data = bl, aes(, x = x, y = bin.id), color = 'white') +
  facet_grid(arm ~ pc.id, scales = 'free', switch = 'y') +
  theme_minimal() +
  geom_point(data = arm.imp, aes(x = abs.value, y = bin.id, color = sign.value)) +
  geom_segment(data = arm.imp, aes(x = abs.value, xend = 0, y = bin.id, yend = bin.id, color = sign.value), linewidth = 1) +
  scale_color_manual(values = cols, na.value = '#999999') +
  labs(x = 'Log odds\nratio', y = '', color = '') +
  theme(strip.text.x = element_blank(),
        strip.text.y.left = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "white", fill = NA, linewidth = 0.2),
        panel.grid.minor.y  = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title.x = element_text(size = 7),
        panel.spacing=unit(0.1, "lines"),
        legend.position = 'none') 


tcga.fig.data <- readRDS(url("https://github.com/cancer-genomics/reproduce_lucas_wflow/raw/master/data/TCGA_Lung/fig2c_p2_data.rds", method="libcurl")) 
tcga.fig.data <- tcga.fig.data %>% mutate(disease = ifelse(disease == "LUAD", "TCGA\nLUAD\n(n = 518)", "TCGA\nLUSC\n(n = 501)"))
tcga.fig.data$disease.new <- factor(tcga.fig.data$disease, levels = c("TCGA\nLUAD\n(n = 518)", "TCGA\nLUSC\n(n = 501)"), 
                        labels = c("TCGA<br>LUAD<br>(*n* = 518)", "TCGA<br>LUSC<br>(*n* = 501)"))

library(data.table)
setDT(tcga.fig.data)
tcga.fig.data[,bin:=as.factor(rev(bin))][]
tcga.plot <- ggplot(tcga.fig.data, aes(x=value, y=bin, color=change, fill=change)) +
  geom_col() +
  facet_grid(arm ~ disease.new, scales = 'free_y', switch = 'y')+
  scale_x_continuous(labels = c(-0.5, 0, 0.5), breaks = c(-0.5, 0, 0.5), minor_breaks = c(-0.2, 0.2)) +
  scale_fill_brewer(palette = 'Set1')+
  scale_color_brewer(palette = 'Set1')+
  theme_minimal() +
  xlab(label = "Proportion of\ncases with CNV") + 
  ylab(label = "") + 
  theme(strip.text.y.left = element_blank(),
        #strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6),
        strip.text.x = ggtext::element_markdown(angle = 0, hjust = 0.5, vjust = 0.5, size = 6),
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 6),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid =  element_blank(),
        legend.title = element_blank(),
        panel.spacing=unit(0.1, "lines"),
        axis.title.x = element_text(size = 7)) + 
  theme(legend.position = "none"
        #legend.text = element_text(size = 6), 
        #legend.key.size = unit(0.4, "cm")
  )+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 

#########################
#then plot actual z-scores from classifier training cohort
#########################

midpoint <- locs %>% 
  group_by(arm) %>% summarize(min = min(bin), max = max(bin)) %>% 
  ungroup() %>% 
  mutate(midpoint = (min + max)/2) %>% 
  select(arm, midpoint) %>% 
  mutate(arm2 = factor(arm, levels = c("1p", "1q", "2p", "2q", "3p", "3q", 
                                                 "4p", "4q", "5p", "5q", "6p", "6q", 
                                                 "7p", "7q", "8p", "8q", "9p", "9q", 
                                                 "10p", "10q", "11p", "11q", "12p", "12q",
                                                 "13q",  "14q", "15q", 
                                                 "16p", "16q", "17p", "17q", "18p", "18q", 
                                                 "19p", "19q", "20p", "20q", "21q", 
                                                 "22q")))


#read in clinical data 
download.file(url = "https://github.com/cancer-genomics/reproduce_lucas_wflow/raw/master/data/LUCAS_metadata.xlsx", destfile = "LUCAS_metadata.xlsx", mode = "wb", quiet = TRUE)
lucas.meta <- readxl::read_excel("LUCAS_metadata.xlsx")

zscores <- features %>% select(id, contains("z")) %>% 
  gather(key = "arm2", value = "zscore", -id) %>% 
  mutate(arm = str_replace(arm2, "zscore_", "")) %>% 
  inner_join(lucas.meta %>% select(id, hist = `Histological diagnosis`)) %>% 
  mutate(hist2 = ifelse(hist %in% c("No baseline cancer", "Benign"), "Noncancer", hist)) %>%
  filter(hist2 %in% c("Noncancer", "Adenocarcinoma", "Squamous"))

#annotation on z-scores for plotting
zs.to.plot <- zscores %>% 
  group_by(arm) %>% 
  mutate(armmax = max(abs(zscore))) %>%
  ungroup() %>% 
  mutate(div = abs(zscore/armmax)) %>%
  mutate(transp = ifelse(div > 0.6, div, 0.6)) %>%
  mutate(hist2 = ifelse(hist2 == "Squamous", "Squamous cell\ncarcinoma", hist2)) %>%
  group_by(hist2, arm) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  mutate(hist.label = paste0(hist2, "\n(n = ", count, ")")) %>%
  select(id, arm, zscore, armmax, transp, id, hist2, count, hist.label) %>%
  mutate(hist.label.f = factor(hist.label, levels = c("Adenocarcinoma\n(n = 62)", "Squamous cell\ncarcinoma\n(n = 29)", "Noncancer\n(n = 158)"), 
                                          labels = c("Adenocarcinoma<br>(*n* = 62)", "Squamous cell<br>carcinoma<br>(*n* = 29)", "Noncancer<br>(*n* = 158)"))) %>%
  mutate(gain.loss = ifelse(zscore > 0, "Gain", "Loss")) %>% 
  inner_join(midpoint)
  


zs.plot <- ggplot() + 
  geom_point(data = zs.to.plot, aes(x=midpoint, y=zscore, color=gain.loss, alpha=transp)) +
  facet_grid(arm2 ~ hist.label.f, scales = 'free_y', switch = 'y') +
  coord_flip() +
  theme_minimal() + 
  theme(strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.ticks = element_blank(),
        panel.grid =  element_blank(),
        plot.title = element_blank(),
        legend.title = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_text(size = 7),
        panel.spacing=unit(0, "lines"),
        strip.text.x = ggtext::element_markdown(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, 'cm'), 
        legend.margin = margin(0,0,0,-10), 
        legend.box.margin = margin(0,0,0,-10)) + 
  scale_y_continuous( breaks = c(-200, 0, 200), minor_breaks = c(-200, 0,200)) + 
  scale_alpha_identity() +
  scale_fill_brewer(palette = 'Set1')+
  scale_color_brewer(palette = 'Set1')+
  ylab("Z-score")+ 
  theme(plot.margin = unit(c(t = 0, b = 0, r = 0, l = 0), "cm")) + 
  guides(colour=guide_legend(override.aes=list(shape=15, size = 5)))



#let's try to bind everything together
#heatmap
#then side annotation
#then TCGA

library(cowplot)
bottom <- plot_grid(tcga.plot, zs.plot, heatmap, right.annot,  align = "h", nrow = 1, rel_widths = c(2.5,4,4.25,1.75))
ggsave(bottom, filename = "bottom.panels.pdf", height = 140, width = 200, units = "mm")
ggsave(top.annot, filename = "top.annotation.pdf", height = 30, width = 67, units = "mm")



