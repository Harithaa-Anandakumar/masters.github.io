library(ggplot2)
library(tidyverse)
library(viridis)
library(patchwork)
#read in data after single cell clustering and group naming to verify that the clusters are enriched for terms similar to that of the paper
#Figure 3.2 

day15group1 <- read.table("data/grp1day30.txt", sep="\t", header = T)

day15group2 <- read.table("data/grp2day30.txt", sep="\t", header = T)

split2 <- day15group2 %>%
  separate(Term, c("GOnum","Process"),"~")
split <- day15group1 %>%
  separate(Term, c("GOnum","Process"),"~")

#GO terms
monkey <- c("sarcomere","system development",
            "cardiovascular system development",
            "cell differentiation","extracellular matrix",
            "extracellular matrix organization",
            "focal adhesion","cell migration",
            "cardiac muscle tissue development",
            "muscle contraction","cardiac muscle contraction",
            "heart development")

daythirty <- c("muscle contraction",
               "sarcomere",
               "muscle tissue development",
               "heart development",
               "system development",
               "animal organ development",
               "cell differentiation",
               "organ morphogenesis",
               "extracellular matrix",
               "cell differentiation",
               "cardiac muscle contraction",
               "sarcomere organization",
               "regulation of cardiac conduction",
               "cardiac myofibril assembly",
               "cardiac muscle tissue morphogenesis",
               "myofibril assembly",
               "cardiac muscle tissue development",
               "adult heart development","extracellular matrix organization",
               "ossification","cell-matrix adhesion","regulation of cell migration")


split <- split %>%
  select(Process, Fold.Enrichment) %>%
  filter(Process %in% daythirty)
split2 <- split2 %>%
  select(Process, Fold.Enrichment) %>%
  filter(Process %in% daythirty)




split$Process <- as.factor(split$Process)
split$group <- rep("one",10)
split2$group <- rep("two",4)


split <- melt(split)
split2 <- melt(split2)
bound <- rbind(split,split2)

thirtyLabs <- c("definitive CM","non-contractile")
fifteenLabs <- c("committed CM","non-contractile")
thirty <- ggplot(bound, aes(x=Process, y=group, size=value, color=group)) +
  geom_point(alpha=0.9) +
  scale_size(range = c(3, 7))+
  xlab("") +
  ylab("")+
  ggtitle("Day 30") +
  guides(colour="none")+
  scale_color_viridis_d()+
  coord_flip()+
  scale_y_discrete(labels=thirtyLabs)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(colour="grey9",angle=45, hjust=1)) +
  guides(color=FALSE, size=FALSE)

  
  
fifteen <- fifteen +
  xlab("") +
  ylab("")+
  ggtitle("Day 15") +
  guides(colour="none")+
  scale_size(range = c(3, 7))+
  scale_color_viridis_d()+
  coord_flip()+
  theme_bw()+
  scale_y_discrete(labels=fifteenLabs)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(colour="grey9",angle=45, hjust=1)) +
  guides(size = guide_legend("Fold Enrichment", override.aes=list(colour="grey70")))



(fifteen | thirty) +
  plot_layout(guides = "collect")
tiff("01_MM.png",units="in", width=8,height=5, res=300, compression = 'lzw') 




