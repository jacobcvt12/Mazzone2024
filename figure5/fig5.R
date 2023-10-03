library(tidyverse)
library(cowplot)

base.size <- 15
font.size <- 3

# figure 5B

data.5b <- tibble(`Screening Scenario`=factor(c("LDCT Alone",
                                                "Blood Test Low Uptake",
                                                "Blood Test High Uptake"),
                                           levels=c("LDCT Alone",
                                                    "Blood Test Low Uptake",
                                                    "Blood Test High Uptake")),
               `Number of cancers detected`=c(24548,
                                              64639,
                                              102431))

fig5b <- ggplot(data.5b, aes(x=`Screening Scenario`,
                          y=`Number of cancers detected`)) +
    geom_bar(stat="identity", aes(fill=`Screening Scenario`),
             show.legend=FALSE) +
    geom_text(aes(label=scales::comma(`Number of cancers detected`)), 
              colour="white", vjust=1.5, size=font.size) +
    geom_segment(aes(x=c(1), y=24548, 
                     xend=c(1), yend=120000)) +
    geom_segment(aes(x=c(1), y=100000, 
                     xend=c(2), yend=100000)) +
    geom_segment(aes(x=c(1), y=120000, 
                     xend=c(3), yend=120000)) +
    geom_segment(aes(x=c(2), y=100000, 
                     xend=c(2), yend=64639)) +
    geom_segment(aes(x=c(3), y=120000, 
                     xend=c(3), yend=102431)) +
    geom_text(label="163% increase", x=1.5, y=95000, size=font.size) +
    geom_text(label="317% increase", x=2, y=115000, size=font.size) +
    theme_classic(base_size=base.size) +
    theme(axis.text=element_text(size=8)) +
    scale_y_continuous(labels=scales::comma, n.breaks=9) +
    scale_fill_manual(values=c("grey", "lightblue", "darkblue")) +
    expand_limits(y=c(0, 160000)) +
    labs()#title="Screen Detected Cancers")

# figure 5C
data.5c <- tibble(`Screening Scenario`=factor(c("LDCT Alone",
                                                "Blood Test Low Uptake",
                                                "Blood Test High Uptake",
                                                "LDCT Alone",
                                                "Blood Test Low Uptake",
                                                "Blood Test High Uptake"),
                                           levels=c("LDCT Alone",
                                                    "Blood Test Low Uptake",
                                                    "Blood Test High Uptake")),
               `Stage at detection`=c(rep("% Detected Stage I", 3),
                                      rep("% Detected Stage IV", 3)),
               `Percent cancers detected`=c(0.262, 0.311, 0.362,
                                            0.450, 0.409, 0.363))

fig5c <- ggplot(data.5c, aes(x=`Stage at detection`,
                    y=`Percent cancers detected`)) +
    geom_bar(stat="identity", position=position_dodge(width=0.9),
             aes(fill=`Screening Scenario`)) +
    geom_text(#x=c(0.8, 1, 1.2, 1.8, 2, 2.2),
              aes(label=scales::percent(`Percent cancers detected`),
                  group=`Screening Scenario`), 
              position=position_dodge(width=0.9),
              colour="white", vjust=1.5, size=font.size) +
    theme_classic(base_size=base.size) +
    theme(legend.position="top") +
    theme(axis.text=element_text(size=8)) +
    scale_y_continuous(labels=scales::percent, n.breaks=6) +
    scale_fill_manual(values=c("grey", "lightblue", "darkblue")) +
    expand_limits(y=c(0, 0.5)) +
    labs(#title="Cancers Detected by Stage and Screening Scenario",
         fill="")

# figure 5D
data.5d <- tibble(`Screening Scenario`=factor(c("LDCT Alone",
                                                "Blood Test Low Uptake",
                                                "Blood Test High Uptake"),
                                           levels=c("LDCT Alone",
                                                    "Blood Test Low Uptake",
                                                    "Blood Test High Uptake")),
               `Lung Cancer Deaths Prevented`=c(4610,
                                                8614,
                                                12866))

fig5d <- ggplot(data.5d, aes(x=`Screening Scenario`,
                             y=`Lung Cancer Deaths Prevented`)) +
    geom_bar(stat="identity", aes(fill=`Screening Scenario`),
             show.legend=FALSE) +
    geom_text(aes(label=scales::comma(`Lung Cancer Deaths Prevented`)), 
              colour="white", vjust=1.5, size=font.size) +
    theme_classic(base_size=base.size) +
    theme(axis.text=element_text(size=8)) +
    scale_y_continuous(labels=scales::comma, n.breaks=3) +
    scale_fill_manual(values=c("grey", "lightblue", "darkblue")) +
    expand_limits(y=c(0, 15000)) +
    labs()#title="Deaths")

# figure 5E
data.5e <- tibble(`Screening Scenario`=factor(c("LDCT Alone",
                                                "Blood Test Low Uptake",
                                                "Blood Test High Uptake"),
                                           levels=c("LDCT Alone",
                                                    "Blood Test Low Uptake",
                                                    "Blood Test High Uptake")),
               `Number needed to scan`=c(202,
                                         164,
                                         156))

fig5e <- ggplot(data.5e, aes(x=`Screening Scenario`,
                          y=`Number needed to scan`)) +
    geom_bar(stat="identity", aes(fill=`Screening Scenario`),
             show.legend=FALSE) +
    geom_text(aes(label=scales::comma(`Number needed to scan`)), 
              colour="white", vjust=1.5, size=font.size) +
    geom_segment(aes(x=c(1), y=202, 
                     xend=c(1), yend=250)) +
    geom_segment(aes(x=c(1), y=220, 
                     xend=c(2), yend=220)) +
    geom_segment(aes(x=c(1), y=250, 
                     xend=c(3), yend=250)) +
    geom_segment(aes(x=c(2), y=220, 
                     xend=c(2), yend=164)) +
    geom_segment(aes(x=c(3), y=250, 
                     xend=c(3), yend=156)) +
    geom_text(label="-38", x=1.5, y=210, size=font.size) +
    geom_text(label="-46", x=2, y=240, size=font.size) +
    theme_classic(base_size=base.size) +
    theme(axis.text=element_text(size=8)) +
    scale_y_continuous(labels=scales::comma, n.breaks=6) +
    scale_fill_manual(values=c("grey", "lightblue", "darkblue")) +
    expand_limits(y=c(0, 250)) +
    labs()#title="Number Needed to Scan")

# assemble
pdf("fig5.pdf", width=10)
plot_grid(fig5b, fig5c, 
          fig5d, fig5e,
          ncol=2)
dev.off()


