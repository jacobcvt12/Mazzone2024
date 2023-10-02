library(aws.s3)
library(tidyverse)
library(cowplot)
library(grid)
library(forcats)

status <- s3read_using(FUN=read_csv,
                       bucket="dlcst-detection-training-crispi",
                       object="features/ldt_v1_training_ids.csv") %>%
  filter(set %in% c("TRAINING", "CV125"))

feats <- s3read_using(FUN=read_csv,
                      bucket="dlcst-detection-training-crispi",
                      object="features/ldt_v1_features_long.csv.gz") %>%
  filter(feature=="centeredslratio") %>%
  inner_join(status) 

nodules <- s3read_using(FUN=read_csv,
                        bucket="delfi-unblinded",
                        object="ldt_analysis/training/clinical_data/training_v1/adsl_delfi-l101_ldt_train_analysis_30MAY2023.csv")

arms <- s3read_using(FUN=read_tsv,
                     bucket="delfi-central-data",
                     object="bins_100Kb_to_5Mb.tsv.gz") %>%
  distinct(armlevel, bin)

arm <- arms %>% group_by(armlevel) %>%
  summarize(n=n(), .groups="drop") %>%
  mutate(arm = as.character(armlevel))

res <- feats %>%
  inner_join(status) %>%
  filter(set=="TRAINING") %>%
  inner_join(nodules) %>%
  # filter(TESTFL) %>%
  mutate(group=case_when(status == "Cancer" ~ "Lung cancer",
                         TRUE ~ "Non-cancer"),
         bin=as.numeric(spread_var),
         group=factor(group, levels=c("Non-cancer",
                                      "Lung cancer"))) %>%
  inner_join(arms) %>%
  mutate(armlevel=fct_reorder(armlevel, bin)) 


text.dat <- res %>%
  group_by(group) %>%
  summarise(N=length(unique(id))) %>%
  ungroup %>%
  mutate(label=paste0(group, " (n=", N, ")"),
         x=-Inf, y=Inf)

x.align <- 0.15

tgNC  <- textGrob(text.dat$label[1],
                  x=unit(x.align, "npc"),
                  y=unit(0.93, "npc"), 
                  gp=gpar(cex=0.5, hjust=0))


tgLC  <- textGrob(text.dat$label[2],
                  x=unit(x.align, "npc"),
                  y=unit(0.48, "npc"), 
                  gp=gpar(cex=0.5, hjust=0))


fig <- res %>%
  # filter(group != "Lung cancer") %>%
  ggplot(aes(x=bin, y=value, group=id)) +
  geom_line(color="grey", linewidth=0.1, alpha=0.8) +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  facet_grid(group~armlevel, 
             scales="free_x", 
             space="free_x",
             switch="x") +
  labs(y="Fragmentation Profile", x="") +
  theme_classic(base_size = 6.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.spacing.y = unit(1, "cm", data=NULL),
        plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")) 


gt_fig <- ggplot_gtable(ggplot_build(fig))

widths <- gt_fig$widths
widths_numeric <- as.numeric(widths)
which_facets <- str_detect(as.character(gt_fig$widths), "null")
which_small <- widths_numeric < 28.6
multipliers <- 28.6/widths_numeric
gt_fig$widths[which_facets & which_small] <- gt_fig$widths[which_facets & which_small] * multipliers[which_facets & which_small]

fig <- ggdraw(gt_fig) + draw_grob(tgNC) + draw_grob(tgLC)

pdf("~/Delfi/ForJacob/mathios_fig_a_ct.pdf", height=3.25, width=7.08661)
fig 
dev.off()