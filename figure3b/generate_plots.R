
library(tidyverse)
library(recipes)
library(GenomicRanges)
library(reshape2)
library(aws.s3)
library(data.table)

########################################
# Feature heatmap + importance barplots
########################################

#load prepped recipe
prepped.recipe<- readRDS("prepped_recipe.RDS")

coefficients <- readRDS("coefficient_posterior.csv")

# PC heatmap-----------
# extract the PC loadings
loadings <- prepped.recipe %>% tidy(id = "pca_AeYA4")
loadings <- data.frame(loadings)
loadings <- dcast(loadings, terms ~ component, value.var = 'value')
loadings$bin.id <- as.numeric(gsub('centeredslratio_', '', loadings$terms))
loadings <- loadings[with(loadings, order(bin.id)),]

# identify # PCs that make 95% of total variance
n.pcs.retained <- (prepped.recipe %>% tidy(id = "pca_AeYA4", type = "variance") %>% 
                     filter(terms == "percent variance") %>% mutate(total_variance = cumsum(value)) %>% 
                     filter(total_variance > 95) %>% pull(component))[1]

loadings <- loadings[,c('terms', 'bin.id', sapply(seq(n.pcs.retained), function(x) paste0('PC', x)))]


# now work on visualization of the loadings
pd <- loadings[,! colnames(loadings) %in% c('terms')]
pd <- melt(pd, id.vars = c('bin.id'))
pd$pc.id <- as.numeric(gsub('PC', '', pd$variable))

# annotate pd with chromosome arm and position
bins5mb <- read.delim("hg19_5Mb_bin_coordinates.txt")
names(bins5mb) <- c("chr", "start", "end", "range", "bin", "arm")
locs <- bins5mb
locs$pos <- apply(locs[,c('start', 'end')], 1, mean)
pd <- merge(pd, locs[,c('chr','pos','arm', 'bin')],
            by.x = 'bin.id', by.y = 'bin', all.x = TRUE)
pd$arm <- factor(pd$arm, levels = unique(locs$arm))


# Rank most important features------------
sc <- coefficients %>% 
  gather(key = "feature", value = "coefficient") %>% 
  group_by(feature) %>% summarize(final.coeff = mean(coefficient))
sc$abs.value = abs(sc$final.coeff)
sc$sign.value = factor(sign(sc$final.coeff), levels = c(-1, 1))
sc$feature.type <- sapply(sc$feature, function(x) strsplit(as.character(x), split = '_')[[1]][1])
sc$feature.type <- factor(sc$feature.type, levels = c('zscore', 'centeredslratio', 'intercept', 'mixturemu', 'mixturesigma', 'mixturetheta', 'mtDNA', 'mxabz'))
sc <- sc[with(sc, order(-abs.value, feature.type)),]
sc$feature <- factor(sc$feature, levels = sc$feature)
sc <- subset(sc, final.coeff != 0)

# determine feature order for the heatmap figure
f <- levels(sc$feature)
pd$pc.id.new <- paste0("PC", pd$pc.id)
pd$pc.id.new <- factor(as.character(pd$pc.id.new),
                   levels = paste0("PC", as.character(as.numeric(gsub('centeredslratio_pc_', '', f[grepl('centeredslratio', f)])))))
pd$pc.id <- factor(as.character(pd$pc.id),
                   levels = as.character(as.numeric(gsub('centeredslratio_pc_', '', f[grepl('centeredslratio', f)]))))

sc$feature <- gsub('zscore_z', 'Z ', gsub('centeredslratio_pc_', 'PC ', sc$feature))
sc$feature <- gsub('PC 0', 'PC', sc$feature)
sc$feature <- gsub('Z 0', 'Z ', sc$feature)
sc$feature <- factor(sc$feature , levels = sc$feature)

#---------------- now make plots -----------------#
# visualize feature by PC heatmap
heatmap <- ggplot(pd, aes(x = pc.id, y = bin.id)) +
  facet_grid(arm ~ pc.id.new, scales = 'free', switch = 'y')+
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
# top annotation -> PC ratios importance
ta <- subset(sc, grepl('PC', feature))
ta$pc.id <- gsub('PC', '', gsub(' ','',ta$feature))
ta$sign.value = factor(ta$sign.value, levels = c(-1, 1))
ta <- ta[,c('pc.id', 'abs.value', 'sign.value')]
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
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3), minor_breaks = c(0, 0.1, 0.2, 0.3)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.text.y=element_blank(),
        strip.text.x = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 0.5),
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
arm.imp <- bl %>% group_by(arm) %>% summarize(bin.id = mean(bin.id), x = mean(x))

arm.lookup <- data.frame(arm= unique(bl$arm), arm.num = seq(39))
arm.imp <- arm.imp %>% inner_join(arm.lookup)

arm.coef <- subset(sc, grepl('Z', feature))[,c('feature', 'abs.value', 'sign.value')]
arm.coef$arm.num <- gsub('Z ', '', arm.coef$feature)
arm.coef$feature <- NULL

arm.imp <- merge(arm.imp, arm.coef[,c('arm.num', 'abs.value', 'sign.value')], by = 'arm.num', all.x = TRUE)

bl$pc.id <- 1
right.annot <- ggplot() +
  geom_point(data = bl, aes(, x = x, y = bin.id), color = 'white') +
  facet_grid(arm ~ pc.id, scales = 'free', switch = 'y') +
  theme_minimal() +
  geom_point(data = arm.imp, aes(x = abs.value, y = bin.id, color = sign.value)) +
  geom_segment(data = arm.imp, aes(x = abs.value, xend = 0, y = bin.id, yend = bin.id, color = sign.value), linewidth = 1) +
  scale_color_manual(values = cols, na.value = '#999999') +
  scale_x_continuous(breaks = c(0, 0.5, 1), minor_breaks = waiver()) +
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


########################################
# Arm-level changes plot with TCGA
########################################
#do the tcga reference data plot
#read in reference bins
ref.bins <- GRanges(read.delim("kb100_genomic_coord_ref.txt") %>% 
  filter(!is.na(bin5mb)) %>%
  group_by(bin5mb) %>% 
  summarize(chr = unique(chr), 
            start = min(start), 
            end = max(end), 
            arm = unique(arm)) %>% 
    select(chr, start, end, arm, bin = bin5mb))

#--------------------------------------------------------#
#initialize functions for TCGA data analyses
get_arm_ranges <- function(assembly){
  
  mySession <- rtracklayer::browserSession()
  GenomeInfoDb::genome(mySession) <- assembly
  
  gaps <- rtracklayer::getTable(rtracklayer::ucscTableQuery(mySession, track="cytoBand"))
  gaps$arm <- substr(gaps$name, 1, 1)
  gaps$id <- apply(gaps[,c('chrom', 'arm')], 1, function(x) gsub('chr', '', paste0(x[1], x[2])))
  
  arm.ranges <- plyr::ddply(gaps, plyr::.(id), plyr::summarize, chrom = unique(chrom), start = min(chromStart), end = max(chromEnd))
  clean.arm.ranges <- subset(arm.ranges,! id %in% c('13p', '14p', '15p', '21p', '22p', 'Xp', 'Xq', 'Yp', 'Yq'))
  clean.arm.ranges <- clean.arm.ranges[,c('chrom', 'start', 'end', 'id')]
  colnames(clean.arm.ranges) <- c('seqnames', 'start', 'end', 'arm')
  return(sort(GRanges(clean.arm.ranges)))
}
find_covered_intervals <- function(intervals, data, min.coverage){
  o <- findOverlaps(intervals, data)
  ow <- width(pintersect(intervals[queryHits(o)], data[subjectHits(o)]))
  base <- intervals[queryHits(o)]
  base$covered.fraction <- ow / width(intervals[queryHits(o)])
  total.feature.coverage <- unlist(lapply(split(base, f = base$bin), FUN = function(x) sum(x$covered.fraction)))
  pos <- intervals[intervals$bin %in% names(which(total.feature.coverage >= min.coverage))]
  return(pos)
}
find_overlapping_genes <- function(ranges, genes, within = 0){
  if (within == 1){
    hits <- unique(as.character(subsetByOverlaps(genes, ranges, type = 'within')$gene_name))
  }
  if (within == 0){
    hits <- unique(as.character(subsetByOverlaps(genes, ranges, type = 'any')$gene_name))
  }
  return(paste(hits, collapse = ','))
}
vectorize <- function(bins, intervals){
  pos.ids <- find_covered_intervals(bins, intervals, 0.9)$bin
  out <- rep(0, length(bins))
  out[bins$bin %in% pos.ids] <- 1
  names(out) <- bins$bin
  return(out)
}
#--------------------------------------------------------#
# determine the cn profile of LUAD tumors
# threshold to call level changes is based on Davoli et al., Science, supp table S1
# TCGA data retrieved from https://github.com/cancer-genomics/reproduce_lucas_wflow/tree/master/data/TCGA_Lung

luad.thresh <- 0.175

luad <- read.delim("tcga_LUAD.tsv")
luad$gain <- ifelse(luad$Segment_Mean > luad.thresh, 1, 0)
luad$loss <- ifelse(luad$Segment_Mean < (-1) * luad.thresh, 1, 0)

luad <- split(luad, f = luad$Sample)
luad <- lapply(luad, function(x) {y = x
y$Sample = NULL
y <- subset(y, Chromosome %in% seq(1,22))
y$Chromosome <- sapply(y$Chromosome, function(x) paste0('chr',x))
return(y)})
luad <- lapply(luad, GRanges)

luad.losses <- lapply(luad, function(x) vectorize(ref.bins, subset(x, loss == 1)))
luad.losses <- do.call(rbind, luad.losses)
aliquot <- sapply(rownames(luad.losses), function(x) strsplit(as.character(x), split = '-')[[1]][4])
luad.loss <- apply(luad.losses[aliquot %in% c('01A', '01B', '02A'),], 2, mean)

# 518 luad
luad.gains <- lapply(luad, function(x) vectorize(ref.bins, subset(x, gain == 1)))
luad.gains <- do.call(rbind, luad.gains)
aliquot <- sapply(rownames(luad.gains), function(x) strsplit(as.character(x), split = '-')[[1]][4])
luad.gain <- apply(luad.gains[aliquot %in% c('01A', '01B', '02A'),], 2, mean)
#--------------------------------------------------------#

#and do the same for the lusc tumors
lusc.thresh <- 0.225

lusc <- read.delim("tcga_LUSC.tsv")
lusc$gain <- ifelse(lusc$Segment_Mean > lusc.thresh, 1, 0)
lusc$loss <- ifelse(lusc$Segment_Mean < (-1) * lusc.thresh, 1, 0)

lusc <- split(lusc, f = lusc$Sample)
lusc <- lapply(lusc, function(x) {y = x
y$Sample = NULL
y <- subset(y, Chromosome %in% seq(1,22))
y$Chromosome <- sapply(y$Chromosome, function(x) paste0('chr',x))
return(y)})
lusc <- lapply(lusc, GRanges)

# 501 LUSC
lusc.losses <- lapply(lusc, function(x) vectorize(ref.bins, subset(x, loss == 1)))
lusc.losses <- do.call(rbind, lusc.losses)
aliquot <- sapply(rownames(lusc.losses), function(x) strsplit(as.character(x), split = '-')[[1]][4])
lusc.loss <- apply(lusc.losses[aliquot %in% c('01A', '01B'),], 2, mean)

lusc.gains <- lapply(lusc, function(x) vectorize(ref.bins, subset(x, gain == 1)))
lusc.gains <- do.call(rbind, lusc.gains)
aliquot <- sapply(rownames(lusc.gains), function(x) strsplit(as.character(x), split = '-')[[1]][4])
lusc.gain <- apply(lusc.gains[aliquot %in% c('01A', '01B'),], 2, mean)

ref.bins$`LUAD Loss` = luad.loss
ref.bins$`LUAD Gain` = luad.gain
ref.bins$`LUSC Loss` = lusc.loss
ref.bins$`LUSC Gain` = lusc.gain

ref.bins <- data.frame(ref.bins)
ref.bins$pos = apply(ref.bins[,c('start', 'end')], 1, mean)
ref.bins[,c('start', 'end', 'width', 'strand')] = NULL
b = melt(ref.bins, id.vars = c('seqnames', 'pos', 'arm' ,'bin'))
b$arm <- factor(b$arm, levels = unique(ref.bins$arm))
b$disease <- sapply(b$variable, function(x) strsplit(as.character(x), split = '[.]')[[1]][1])
b$change <- sapply(b$variable, function(x) strsplit(as.character(x), split = '[.]')[[1]][2])
b <- b %>% mutate(disease = ifelse(disease == "LUAD", "TCGA\nLUAD\n(n = 518)", "TCGA\nLUSC\n(n = 501)"))
b$disease.new <- factor(b$disease, levels = c("TCGA\nLUAD\n(n = 518)", "TCGA\nLUSC\n(n = 501)"), 
                                 labels = c("TCGA<br>LUAD<br>(*n* = 518)", "TCGA<br>LUSC<br>(*n* = 501)"))
b[which(b$change == 'Loss'), 'value'] = b[which(b$change == 'Loss'), 'value'] * (-1)
#--------------------------------------------------------#

tcga.fig.data <- b

setDT(tcga.fig.data)
tcga.plot <- ggplot(tcga.fig.data, aes(x=value, y=bin, color=change, fill=change)) +
  facet_grid(arm ~ disease.new, scales = 'free_y', switch = 'y')+
  geom_col() +
  scale_x_continuous(labels = c(-0.5, 0, 0.5), breaks = c(-0.5, 0, 0.5), minor_breaks = c(-0.2, 0.2)) +
  scale_fill_brewer(palette = 'Set1')+
  scale_color_brewer(palette = 'Set1')+
  theme_minimal() +
  xlab(label = "Proportion of\ncases with CNV") + 
  ylab(label = "") + 
  theme(strip.text.y.left = element_blank(),
    strip.text.x = ggtext::element_markdown(angle = 0, hjust = 0.5, vjust = 0.5, size = 6),
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 6),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid =  element_blank(),
        legend.title = element_blank(),
        panel.spacing=unit(0.1, "lines"),
        axis.title.x = element_text(size = 7)) + 
  theme(legend.position = "none")+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 


#########################
#then plot actual z-scores from classifier training cohort
#########################
#read in sample IDs and clinical info
sample_ids <- read.csv("training_manifest.csv")

#gather z-scores
zscores <- readRDS("training_z_scores.rds")

#map chromosome arm to bin
arm.map <- read.delim("zscores_bins.txt") %>% select(arm = z.arm.pad, armlevel = chr.arm) %>% unique() %>% 
  mutate(armlevel = factor(armlevel, levels = c("1p", "1q", "2p", "2q", "3p", "3q", 
                                                 "4p", "4q", "5p", "5q", "6p", "6q", 
                                                 "7p", "7q", "8p", "8q", "9p", "9q", 
                                                 "10p", "10q", "11p", "11q", "12p", "12q",
                                                  "13q",  "14q", "15q", 
                                                 "16p", "16q", "17p", "17q", "18p", "18q", 
                                                 "19p", "19q", "20p", "20q", "21q", 
                                                  "22q")))

#annotation on z-scores for plotting
zs.to.plot <- zscores %>% 
  group_by(arm) %>% 
  mutate(armmax = max(abs(z))) %>%
  ungroup() %>% 
  mutate(div = abs(z/armmax)) %>%
  mutate(transp = ifelse(div > 0.6, div, 0.6)) %>%
  mutate(id = sample_id) %>% 
  inner_join(sample_ids) %>% 
  mutate(hist = ifelse(LCNCAT %in% c("Not Confirmed Lung Cancer"), "Noncancer", LCNHIST)) %>% 
  mutate(hist = ifelse(hist == "Noncancer", "Noncancer", 
                       ifelse(LCNHIST %in% c("Adenocarcinoma", "Squamous cell carcinoma"), LCNHIST, "other"))) %>% 
  mutate(hist = ifelse(hist == "Squamous cell carcinoma", "Squamous cell\ncarcinoma", hist)) %>%
  group_by(LCNHIST, arm) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  mutate(hist.label = paste0(hist, "\n(n = ", count, ")")) %>%
  select(sample_id, arm, z, armmax, transp, id, hist, count, hist.label) %>%
  mutate(hist.label.f = factor(hist.label, levels = c("Adenocarcinoma\n(n = 89)", "Squamous cell\ncarcinoma\n(n = 57)", "Noncancer\n(n = 395)"), 
                                          labels = c("Adenocarcinoma<br>(*n* = 89)", "Squamous cell<br>carcinoma<br>(*n* = 57)", "Noncancer<br>(*n* = 395)"))) %>%
  mutate(gain.loss = ifelse(z > 0, "Gain", "Loss")) %>% 
  inner_join(arm.map) %>% filter(!is.na(hist.label.f))
  

zs.plot <- ggplot() + 
  geom_point(data = zs.to.plot, aes(x=armlevel, y=z, color=gain.loss, alpha=transp)) +
  facet_grid(armlevel ~ hist.label.f, scales = 'free_y', switch = 'y') +
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
  scale_y_continuous(limits = c(-150,150)) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, 'cm'), 
        legend.margin = margin(0,0,0,-10), 
        legend.box.margin = margin(0,0,0,-10)) + 
  scale_alpha_identity() +
  scale_fill_brewer(palette = 'Set1')+
  scale_color_brewer(palette = 'Set1')+
  ylab("Z-score")+ 
  theme(plot.margin = unit(c(t = 0, b = 0, r = 0, l = 0), "cm")) + 
  guides(colour=guide_legend(override.aes=list(shape=15, size = 5)))



#let's try to bind everything together
library(cowplot)

bottom <- plot_grid(tcga.plot, zs.plot, heatmap, right.annot,  align = "h", nrow = 1, rel_widths = c(2.5,4,4.25,1.75))
ggsave(bottom, filename = "bottom.panels.pdf", height = 140, width = 200, units = "mm")
ggsave(top.annot, filename = "top.annotation.pdf", height = 30, width = 67, units = "mm")

