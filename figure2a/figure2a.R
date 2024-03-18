library(aws.s3)
library(tidyverse)
library(cowplot)
library(grid)
library(forcats)
library(RCurl)

# Load preds to get training data status
load(url("https://raw.githubusercontent.com/cancer-genomics/reproduce_lucas_wflow/master/code/rlucas/data/prediction_lucas.rda"))
status <- select(preds, id, status=type)  

# Calculate ratios for 473 bins
feats <-
  read.csv(text = getURL("https://raw.githubusercontent.com/cancer-genomics/reproduce_lucas_wflow/master/data/lucas_5mbs_delfi473.csv")) %>% 
  filter(!str_detect(arm, "X")) %>%
  mutate(ratio.cor = short.cor/long.cor) %>%
  group_by(id) %>%
  mutate(ratio.centered = scale(ratio.cor, scale=FALSE)[,1]) %>%
  ungroup()

arms <- distinct(feats, bin, arm) %>%
  rename(armlevel=arm)

arm <- arms %>% group_by(armlevel) %>%
  summarize(n=n(), .groups="drop") %>%
  mutate(arm = as.character(armlevel))

# Create data frame to plot
res <- feats %>%
  inner_join(status) %>%
  mutate(group=case_when(status == "cancer" ~ "Lung cancer",
                         TRUE ~ "Non-cancer"),
         group=factor(group, levels=c("Non-cancer",
                                      "Lung cancer"))) %>%
  inner_join(arms) %>%
  mutate(armlevel=fct_reorder(armlevel, bin)) 

# Create data frame for text grobs
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


# Make the figure
fig <- res %>%
  ggplot(aes(x=bin, y=ratio.centered, group=id)) +
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

# Save figure as PDF
pdf("figure2a.pdf", height=3.25, width=7.08661)
fig 
dev.off()