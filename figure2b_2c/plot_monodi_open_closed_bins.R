library(tidyverse)
library(aws.s3)
library(scales)
set.seed(1207)


FindOvls <- function(obj1,obj2){
  
  keepi <- findOverlaps(obj1,obj2)
  freq.matched <- obj1[queryHits(keepi)]
  
  mcols(freq.matched) <- cbind.data.frame(
    mcols(freq.matched),
    mcols(obj2[subjectHits(keepi)]))
  return(as.data.frame(freq.matched))
}

.meanSmoother <- function(x, k=1, iter=2, na.rm=TRUE){
  meanSmoother.internal <- function(x, k=1, na.rm=TRUE){
    n <- length(x)
    y <- rep(NA,n)
    
    window.mean <- function(x, j, k, na.rm=na.rm){
      if (k>=1){
        return(mean(x[(j-(k+1)):(j+k)], na.rm=na.rm))
      } else {
        return(x[j])
      }
    }
    
    for (i in (k+1):(n-k)){
      y[i] <- window.mean(x,i,k, na.rm)
    }
    for (i in 1:k){
      y[i] <- window.mean(x,i,i-1, na.rm)
    }
    for (i in (n-k+1):n){
      y[i] <- window.mean(x,i,n-i,na.rm)
    }
    y
  }
  
  for (i in 1:iter){
    x <- meanSmoother.internal(x, k=k, na.rm=na.rm)
  }
  x
}

#read in fragmentation profiles
scaled.dat <- read_csv("figure2b_2c/Lucas_monodi_ratios_100kb_bins.csv")

#gather supporting reference files
ref.bins <- read.delim("figure2b_2c/kb100_genomic_coord_ref.txt")

#LUSC Reference AB compartments
lusc.data <- read.delim("figure2b_2c/lusc_tumor_compartments_100kb.txt") %>% 
 select(chr = chr, start, end,lusc=eigen)%>%
 merge(ref.bins, by=c("chr", "start", "end")) %>%
 select(c(bin_num, arm, lusc))%>% 
 mutate(lusc=scale(lusc, center=T, scale=T)) %>%
 mutate(lusc=rescale(lusc))

#smooth bin eigenvalues
lusc <- lusc.data %>%
 arrange(bin_num) %>%
group_by(arm) %>%
nest() %>%
mutate(data=map(data,
                ~ mutate(.x,
                         lusc=.meanSmoother(lusc)))) %>%
 unnest("data")%>%
  ungroup()%>%
 select(-c(arm))

lymph.data <- read.delim("figure2b_2c/hic_compartments_100kb_ebv_2014.txt", sep = " ") %>% 
  select(chr = chr, start, end, lymph=eigen)%>%
  merge(ref.bins, by=c("chr", "start", "end")) %>%
  select(c(bin_num, arm, lymph)) %>%
  mutate(lymph=scale(lymph, center=T, scale=T)) %>%
  mutate(lymph=rescale(lymph))

#smooth bin eigenvalues
lymph <- lymph.data  %>%
  arrange(bin_num) %>%
  group_by(arm) %>%
  nest() %>%
  mutate(data=map(data,
                  ~ mutate(.x,
                          lymph=.meanSmoother(lymph)))) %>%
  unnest("data") %>%
  ungroup()%>%
  select(-c(arm))


#lymphoblastoid Reference AB compartments

data <- read.delim("~/Downloads/lusc_tumor_compartments_100kb.txt") %>% 
  select(chr = chr, start, end,lusc=eigen)%>%
  merge(ref.bins, by=c("chr", "start", "end")) %>%
  select(c(bin_num, arm, lusc))
#smooth bin eigenvalues
lusc.dat <- data %>%
  arrange(bin_num) %>%
  group_by(arm) %>%
  nest() %>%
  mutate(data=map(data,
                  ~ mutate(.x,
                           lusc=.meanSmoother(lusc)))) %>%
  unnest("data")%>%
  ungroup() %>%
  select(-c(arm))

data <- read.delim("~/Downloads/hic_compartments_100kb_ebv_2014.txt", sep = " ") %>% 
  select(chr = chr, start, end, lymph=eigen)%>%
  merge(ref.bins, by=c("chr", "start", "end")) %>%
  select(c(bin_num, arm, lymph))

lymph.dat <- data  %>%
  arrange(bin_num) %>%
  group_by(arm) %>%
  nest() %>%
  mutate(data=map(data,
                  ~ mutate(.x,
                           lymph=.meanSmoother(lymph)))) %>%
  unnest("data") %>%
  ungroup() %>%
  select(-c(arm))

#set reference and test non-cancers
ref.noncancers <-read_csv("figure2b_2c/Noncancer_reference_ids.txt")

test.noncancers <- read_csv("figure2b_2c/Noncancer_test_ids.txt") %>%
  select(-c(ichor))

sclc.samps <- read_csv("figure2b_2c/Cancer_ids.txt")

samps <- rbind(sclc.samps, test.noncancers)

#generate reference value of each metric across 10 healthy cases
reference_healthies <- scaled.dat%>%
  filter(id %in% ref.noncancers$id) %>%
  group_by(metric, bin_num) %>%
  summarise(ref_healthies=mean(scaled.smooth))

#generate reference open & closed regions
#open in lusc, closed in lymgh

refs_lymph_closed <- merge(lymph, lusc, by="bin_num") %>%
  mutate(dif=lusc-lymph) %>%
  arrange(dif) %>%
  filter(lymph >= .7 & lusc <=.3) %>%
  mutate(color="lymph_closed_lusc_open")


refs_lymph_open <- merge(lymph, lusc, by="bin_num") %>%
  mutate(dif=lusc-lymph) %>%
  arrange(dif) %>%
  filter(lymph <= .3 & lusc >=.7) %>%
  mutate(color="lymph_open_lusc_closed")


refs_all_closed <- merge(lymph, lusc, by="bin_num") %>%
  mutate(dif=lusc-lymph) %>%
  arrange(dif) %>%
  filter(dif > -.01  & dif < .01) %>%
  filter(lymph > .65) %>%
  mutate(color="all_closed")


refs_all_open <- merge(lymph, lusc, by="bin_num") %>%
  mutate(dif=lusc-lymph) %>%
  arrange(dif) %>%
  filter(dif > -.01  & dif < .01) %>%
  filter(lymph < .65) %>%
  mutate(color="all_open")


refs <- rbind(refs_lymph_closed, refs_lymph_open) %>%
  rbind(refs_all_closed) %>%
  rbind(refs_all_open) 


dat <- scaled.dat %>%
  filter(!id %in% c(ref.noncancers$id)) %>%
  merge(refs, by="bin_num") %>%
  merge(reference_healthies, by=c("bin_num", "metric")) %>%
  mutate(scaled.smooth=scaled.smooth-ref_healthies) %>%
  group_by(id, color, metric) %>%
  summarise(med=mean(scaled.smooth, na.rm=T)) %>%
  mutate(type = case_when(
    id %in% sclc.samps$id ~ "Cancers", 
    id %in% test.noncancers$id ~ "Non-cancer")
  )%>%
  mutate(region=ifelse(color %in% c("all_closed", "lymph_closed_lusc_open"), "closed", "open"))


#generate summarized metric across cancers, noncancers in 100kb bins
merged.smooth.type <- scaled.dat%>%
  merge(samps, by="id") %>%
  filter(metric=="monodi") %>%
  merge(reference_healthies, by=c("bin_num", "metric")) %>%
  mutate(scaled.smooth=scaled.smooth-ref_healthies) %>%
  filter(type=="Cancer") %>%
  dplyr::rename("Cancer"=scaled.smooth) %>%
  select(-c(type)) %>%
  gather(type, scaled.smooth, c(Cancer, ref_healthies)) %>%
  group_by(bin_num, type) %>%
  summarise(monodi=median(scaled.smooth)) %>%
  gather(variable, value, -c(bin_num, type)) %>%
  unite(temp, type, variable) %>%
  spread(temp, value) %>%
  merge(lymph.dat, by="bin_num") %>%
  merge(lusc.dat, by="bin_num") %>%
  merge(ref.bins, by="bin_num") %>%
  filter(arm != "#N/A") %>%
  gather(metric, value, c(Cancer_monodi:lusc)) %>%
  mutate(color=ifelse(value > 0, "gray50", "red4")) %>%
  na.omit()

highlight <- merged.smooth.type %>%
  select(c(bin_num, metric, color)) %>%
  spread(metric, color) %>%
  mutate(highlight=ifelse(Cancer_monodi == lusc & ref_healthies_monodi==lymph & lusc != lymph, "yes", "no")) %>%
  filter(highlight=="yes") %>%
  select(c(bin_num, highlight))

#fragmentation profile of mono/di ratio in 100kb bins
p1 <- merged.smooth.type %>%
  filter(chr == "chr22") %>%
  group_by(metric) %>%
  mutate(value=as.numeric(value)) %>%
  mutate(value=scale(value, center=F, scale=T)) %>%
  mutate(transp=ifelse(bin_num %in% highlight$bin_num, .95, .1)) %>%
  mutate(metric=as.character(metric),
         metric=str_replace_all(metric, "Cancer_monodi", "LUSC cfDNA Component"),
         metric=str_replace_all(metric, "lusc","LUSC Tissue Reference AB Compartments"),
         metric=str_replace_all(metric,"ref_healthies_monodi","Individuals without cancer cfDNA"),
         metric=str_replace_all(metric,"lymph","Lymphoblastoid cell HiC Reference AB Compartments")) %>%
  mutate(metric=as.factor(metric)) %>%
  mutate(transp=ifelse(bin_num %in% highlight$bin_num, .95, .1)) %>%
ggplot(aes(x=as.factor(bin_num), y=value, fill=color, alpha=transp))+geom_bar(stat="identity")+
  facet_wrap(~factor(metric, levels=c("LUSC Tissue Reference AB Compartments", "LUSC cfDNA Component",  "Individuals without cancer cfDNA", "Lymphoblastoid cell HiC Reference AB Compartments")), ncol=1, scales="free")+
  theme_classic()+ 
  scale_fill_manual(values=c("red4", "gray50"),labels=c("Closed","Open"),name="A/B Compartments")+
  scale_alpha_identity()  +
    xlab("") +
    ylab("Fragmentation Profile") +
    xlab("chr22") +
    theme_classic(base_size=7) +
    theme(legend.position="bottom",
          ##axis.title.y=element_blank(),
          strip.background=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
  scale_y_continuous(breaks=seq(-1,1,1))+theme(legend.position = "top")+
  theme(legend.key.size = unit(.25, 'cm'))

#deviation from non-cancer reference in important bins
p2 <- dat %>%
  filter(metric=="monodi") %>%
  mutate(type.x=as.character(type)) %>%
  mutate(type.x=str_replace(type.x, "Cancers", "Squamous Cell\nLung Cancer" )) %>%
  mutate(region=str_replace(region, "closed", "Closed in \nLymphoblastoid\nReference"), region=str_replace(region, "open", "Open in\nLymphoblastoid\nReference")) %>%
  mutate(color=str_replace(color, "all_open", "Shared Chromatin State"), color=str_replace(color, "all_closed", "Shared Chromatin State"), color=str_replace(color, "lymph_closed_lusc_open", "Different Chromatin State"), color=str_replace(color, "lymph_open_lusc_closed", "Different Chromatin State")) %>%
  ggplot(aes(x=region, y=med, color=type.x, shape=type.x))+labs(x="", y="Median Mono/Di Ratio")+geom_hline(yintercept = 0, linetype="dashed", color="darkgray")+facet_wrap(~color, ncol=2, scales="free_x")+labs(color="", fill="", x="", y="Deviation from non-cancer\ncfDNA reference")+
  geom_point(size=.5, position = position_jitterdodge(jitter.width=.05, jitter.height = .01)) +
  geom_boxplot(aes(fill=type.x), outlier.shape = NA, alpha=0.5) + ##figure out how to
  theme(panel.grid=element_blank(),
        ##panel.spacing=unit(0.5, "lines"),
        ##strip.placement="outside",
        strip.background=element_blank(),
        strip.text=element_text(color="transparent"),
        ##strip.text=element_text(size=22),
        axis.ticks.y=element_blank())+
  scale_color_manual(name="", values=c("gray50", "steelblue"))+
  scale_fill_manual(name="", values=c("gray50", "steelblue"))+
  scale_shape_manual(name="", values=c("circle", "triangle"))+
  theme_classic(base_size=7)+
  theme(legend.key.size = unit(.25, 'cm'))+theme(strip.background=element_blank())

library(cowplot)
plot_final <- plot_grid(p1, p2, rel_heights = c(3,1), ncol=1, align = "v",  axis = 'l', label_size = 12, label_x = .05, hjust=.01)
ggsave("figure2b_2c/figure2bc.pdf", plot_final, units="mm", width=180, height = 170)
