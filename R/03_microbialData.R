#################################

library(readxl)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(numform)
library(viridis)
library(hrbrthemes)
library(purrr)
library(tidyr)
library(dplyr)
library(reshape)
library(ggthemes)
library(ggsci)
library(patchwork)

metdat <- read_xlsx("data/basicStats.xlsx")
metdat <- metdat[-4,]
metdat$unmapped <-metdat$`Total Reads` - rowSums(metdat[,-(1:2)])


metdat <- metdat %>%
  mutate(totalMicrobial = Bacterial+Viral) %>%
  mutate(humanPerc = `Reads Mapped Human`/`Total Reads` * 100) %>%
  mutate(microbialPerc = totalMicrobial/`Total Reads` * 100) %>%
  mutate(rRNAperc = rRNA/`Total Reads` * 100) %>%
  mutate(unmappPerc = unmapped/`Total Reads` * 100)
  #mutate(numbered = seq(1:8))

metdat$Sample <- as.factor(metdat$Sample)
metdat$numbered <- as.factor(metdat$numbered)


onlyPer <- metdat %>%
  select(Sample,contains("erc"))

colnames(onlyPer) <- c("Sample","Primary - Human","Microbial genomes","rRNA Human", "Multimapped" )
onlyPer <- as.data.frame(onlyPer)
onlyPer <- melt(onlyPer)
onlyPer$variable <- as.factor(onlyPer$variable)
levels(onlyPer$variable)
onlyPer$variable <- factor(onlyPer$variable, levels = c("Primary - Human", "Multimapped", "rRNA Human", "Microbial genomes"))
onlyPer$Sample <- as.factor(onlyPer$Sample)

breaks = c(60,30,5,0,0.1,0.01,0.001)
mew <- f_comma(c(60,30,5,0,0.1,.01,.001))
lab <- seq(1:7)

ggplot(onlyPer, (aes(x=Sample, y=(value),colour=variable)))+
  geom_point(stat="identity",na.rm=TRUE,size=4)
m <- ggplot(onlyPer, (aes(x=Sample, y=(value),colour=variable)))+
  geom_point(stat="identity",na.rm=TRUE,size=4) +
  #scale_color_tableau(palette = "Summer")+
  labs(y="% of total reads",
       x= "Samples")+
  theme_bw()+
  labs(face="bold")+
  scale_y_log10(breaks = breaks,
                labels = mew) +
  scale_x_discrete(labels = lab,limits = c("SRR1663123_GSM1554465",
                                           "SRR6706796_GSM2991857",
                                           "Sample_r733sCDICM3",
                                           "p556sCM10-3-4",
                                           "p637sDiff6CM",
                                           "p722s3C190604",
                                           "p786sC190924A"))+
  scale_color_manual(values = c("#9d983d","#90728f","#6b6b6b", "#e15759" ))+
 #scale_color_viridis_d(option =  "plasma", name = NULL) +
  # scale_y_log10(breaks = breaks1, labels = function(x) format(x,scientific=FALSE))+
  theme(legend.key=element_blank(), 
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
        legend.title = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, size = 1.2), 
        legend.position = "right")

m
head(metdat)
eh <- metdat %>%
  mutate(
    PerMilHumMap = `Reads Mapped Human`/10^6,
    numBacPerMilHum = Bacterial/PerMilHumMap,
    numVirPerMilHum = Viral/PerMilHumMap)
eh <- as.data.frame(eh)

eh %>%
  select(Sample, contains("num")) %>%
  melt() %>%
  ggplot()+
  scale_x_discrete(label=lab, limits = c("SRR1663123_GSM1554465",
                                         "SRR6706796_GSM2991857",
                                         "Sample_r733sCDICM3",
                                         "p556sCM10-3-4",
                                         "p637sDiff6CM",
                                         "p722s3C190604",
                                         "p786sC190924A"))+
  geom_point(aes(x=Sample, y=value, colour=variable), size=5)+
  labs(y="Viral / Bacterial reads per million human mapped reads",
       x="Samples")+
  theme_bw()+
  scale_color_manual(values = c("#ff9888", "#cecb76" ),
                      labels = c("Bacteria", "Virus"))+
  theme(text=element_text(size=10))+
  theme(legend.key=element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, size = 1.2), 
        legend.position = "right")

tiff("overviewCon",units="in", width=10,height=4, res=300, compression = 'lzw') 

(m|n) + plot_annotation(tag_levels = 'A')

dev.off()

##############################VIRUS###################

 f_files<- list.files("data/decon/virus",
                     pattern = "*summary*", full.names = T)


read_in_feature_counts<- function(file){
  cnt<- read.table(file, sep="\t", header = T)
  cnt <- cnt %>%
    mutate(fileName=file)
  return(cnt)
}


raw_counts <- map(f_files, read_in_feature_counts)

raw_counts_df<- purrr::reduce(raw_counts, union) 

raw_counts_df <- as.data.frame(raw_counts_df)
raw_counts_df$fileName <- gsub("data/decon/virus/unmapped_","\\", raw_counts_df$fileName)
raw_counts_df$fileName <- gsub("_RIBO_UNALIGNED_vs_viruses_subject_summary_all_CT_5.txt","\\", raw_counts_df$fileName)


metdat$Sample %in% raw_counts_df$fileName
raw_counts_df <- raw_counts_df %>%
                  filter(fileName != "iBM76-3-control_2")
head(raw_counts_df)
colnames(raw_counts_df)
spreaded <- tidyr::spread(raw_counts_df, fileName, Count)

spreaded[is.na(spreaded)] <- 0
orderSample <- metdat$Sample

metdat$Sample %in% colnames(spreaded)
rownames(spreaded) <- spreaded$Species.name
spreaded <- spreaded[,-(1)]
spreaded <- spreaded %>%
  select("SRR1663123_GSM1554465_Fetal_heart_RNA_seq_rep1_Homo_sapiens_RNA-Seq",
         "SRR6706796_GSM2991857_45_HM_3_FC2_Homo_sapiens_RNA-Seq",
         "Sample_r733sCDICM3_hs_totalrna_sr_Zimmermann_D_3",
         "p556sCM10-3-4_Zimmermann_S6_L001_R1_001",
         "p637sDiff6CM_Zimmermann_S46_L005_R1_001",
         "p722s3C190604_Tiburcy_S37_L005_R1_001",
         "p786sC190924A_Tiburcy_S2_L001_R1_001")

## ordering column based on the same order as before 

#doesn't work anymore cuz i messed up metdat 
spreaded <- spreaded[, metdat$Sample]

spreaded$Species.name <- rownames(spreaded)
spreaded <- as.data.frame(spreaded)
viralReads <- colSums(spreaded)

meltd <- melt(spreaded)
head(meltd)

lab <- seq(1:7)



# !!! Create a custom color scale -- keeps colours constant across plots of the same factors 
meltd$Species.name <- gsub("NC_[0-9][0-9][0-9][0-9][0-9][0-9].[0-9]","\\", meltd$Species.name)
meltd$Species.name <- as.factor(meltd$Species.name)
myColors <- tableau_color_pal(palette = "Winter")(7)
names(myColors) <- levels(meltd$Species.name)
colScale <- scale_fill_manual(name = "Species.name",values = myColors)

g <- meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity") +
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Number Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  colScale 

##dplyr not in ! 
meltd$Species.name
subdata <- meltd %>%
  filter(!Species.name %in% c(" Proteus phage"," Col phage"))
  
h <- subdata %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity") +
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Number Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  guides(fill=FALSE)+
  colScale 


i <- subdata %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity", position = "fill") +
  scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Percentage of Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  guides(fill=FALSE)+
  colScale 

j <- subdata %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~Species.name)+
  theme_minimal()+
  labs(x="Samples",
       y="Number of Reads Mapped")+
  scale_x_discrete(label=lab) +
  guides(fill=FALSE)+
  theme(panel.grid.major = element_blank(), 
        panel.border = element_rect(colour = "grey", 
                                    fill = NA, 
                                    size = 1.0),
        strip.text.x = element_text(size=8),
        axis.text.y = element_text(size = 7)
  ) +
  colScale 

pat <- (g + h + i) / j

pat[[1]] <- pat[[1]] + plot_layout(tag_level = 'new',
                                   guides = "collect") 
tiff("viralSp.png",units="in", width=10,height=7, res=300, compression = 'lzw') 
pat + plot_annotation(tag_levels = c('A','1'),
                      tag_sep = '.',
                      title = 'Possible Viral Contaminants')

dev.off()

##############################BACTERIA###################

f<- list.files("data/decon/bacteria",
                     pattern = "*summary_ge*", full.names = T)
read_in_feature_counts<- function(file){
  cnt<- read.table(file, sep="\t", header = T)
  cnt <- cnt %>%
    mutate(fileName=file)
  return(cnt)
}

raw_counts <- map(f, read_in_feature_counts)

raw_counts_df<- purrr::reduce(raw_counts, union) 

raw_counts_df <- as.data.frame(raw_counts_df)
raw_counts_df$fileName <- gsub("data/decon/bacteria/unmapped_","\\", raw_counts_df$fileName)
raw_counts_df$fileName <- gsub("_RIBO_UNALIGNED_vs_bacteria_subject_summary_ge_CT_5.txt","\\", raw_counts_df$fileName)

raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum) %>%
  summarize(sum(Count))


raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum) %>%
  arrange(desc(Count)) %>%
  top_n(10) %>%
  summarize(sum(Count))

##87.7 % contributed by first 10 genus
topTen <- raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum) %>%
  arrange(desc(Count)) %>%
  top_n(10) 

metdat$Sample %in% raw_counts_df$fileName

head(raw_counts_df)

raw_counts_df <- raw_counts_df %>%
  filter(Species.name %in% topTen$Species.name)
 
spreaded <- tidyr::spread(raw_counts_df, fileName, Count)
spreaded[is.na(spreaded)] <- 0

metdat$Sample %in% colnames(spreaded)
rownames(spreaded) <- spreaded$Species.name
spreaded <- spreaded[,-1]
## ordering column based on the same order as before 
spreaded <- spreaded %>%
  select("SRR1663123_GSM1554465_Fetal_heart_RNA_seq_rep1_Homo_sapiens_RNA-Seq",
         "SRR6706796_GSM2991857_45_HM_3_FC2_Homo_sapiens_RNA-Seq",
         "Sample_r733sCDICM3_hs_totalrna_sr_Zimmermann_D_3",
         "p556sCM10-3-4_Zimmermann_S6_L001_R1_001",
         "p637sDiff6CM_Zimmermann_S46_L005_R1_001",
         "p722s3C190604_Tiburcy_S37_L005_R1_001",
         "p786sC190924A_Tiburcy_S2_L001_R1_001")

#spreaded <- spreaded[, metdat$Sample]

spreaded$Species.name <- rownames(spreaded)
meltd <- melt(spreaded)
head(meltd)

lab <- seq(1:7)

meltd$Species.name <- as.factor(meltd$Species.name)

# !!! Create a custom color scale -- keeps colours constant across plots of the same factors 

myColors <- tableau_color_pal(palette = "Winter")(10)
names(myColors) <- levels(meltd$Species.name)
colScale <- scale_fill_manual(name = "Species.name",values = myColors)


x <- meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity") +
  #scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Number Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  colScale

y <-  meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity", position = "fill") +
  scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Percentage of Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  colScale 

z <- meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~Species.name)+
  theme_minimal()+
  labs(x="Samples",
       y="Number of Reads Mapped")+
  scale_x_discrete(label=lab) +
  guides(fill=FALSE)+
  theme(panel.grid.major = element_blank(), 
        panel.border = element_rect(colour = "grey", 
                                    fill = NA, 
                                    size = 1.0),
        strip.text.x = element_text(size=8),
        axis.text.y = element_text(size = 7)
  ) +
  colScale 

pat <- ((x + y)/ z) 

pat[[1]] <- pat[[1]] + plot_layout(tag_level = 'new',
                                   guides = "collect") 
tiff("bacGenus1.tiff",units="in", width=10,height=7, res=300, compression = 'lzw') 
pat + plot_annotation(tag_levels = c('A','1'),
                      tag_sep = '.',
                      title = 'Possible Bacterial Contaminants at Genus Level')

dev.off()


########## species level ##### 

f<- list.files("data/decon/bacteria",
               pattern = "*summary_sp*", full.names = T)
read_in_feature_counts<- function(file){
  cnt<- read.table(file, sep="\t", header = T)
  cnt <- cnt %>%
    mutate(fileName=file)
  return(cnt)
}

raw_counts <- map(f, read_in_feature_counts)

raw_counts_df<- purrr::reduce(raw_counts, union) 


raw_counts_df <- as.data.frame(raw_counts_df)
raw_counts_df$fileName <- gsub("data/decon/bacteria/unmapped_","\\", raw_counts_df$fileName)
raw_counts_df$fileName <- gsub("_RIBO_UNALIGNED_vs_bacteria_subject_summary_sp_CT_5.txt","\\", raw_counts_df$fileName)
dplyr::filter(raw_counts_df, grepl('P', Species.name))

head(raw_counts_df)

raw_counts_df[is.na(raw_counts_df)] <- 0
raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum)  %>%
  summarize(sum(Count))

raw_counts_df %>%
  filter(fileName == "p722s3C190604_Tiburcy_S37_L005_R1_001") %>%
  arrange(desc(Count)) %>%
  select(-fileName)%>%
  top_n(2)%>%
  summarise(sum(Count))


raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum) %>%
  arrange(desc(Count)) %>%
  top_n(20) %>%
  summarize(sum(Count))

##57.7 % contributed by first 10 species
topTen <- raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum) %>%
  arrange(desc(Count)) %>%
  top_n(10) 

metdat$Sample %in% raw_counts_df$fileName

head(raw_counts_df)

raw_counts_df <- raw_counts_df %>%
  filter(Species.name %in% topTen$Species.name)

spreaded <- tidyr::spread(raw_counts_df, fileName, Count)
spreaded[is.na(spreaded)] <- 0

metdat$Sample %in% colnames(spreaded)
rownames(spreaded) <- spreaded$Species.name

spreaded <- spreaded[,-1]
## ordering column based on the same order as before 
spreaded <- spreaded %>%
  select("SRR1663123_GSM1554465_Fetal_heart_RNA_seq_rep1_Homo_sapiens_RNA-Seq",
         "SRR6706796_GSM2991857_45_HM_3_FC2_Homo_sapiens_RNA-Seq",
         "Sample_r733sCDICM3_hs_totalrna_sr_Zimmermann_D_3",
         "p556sCM10-3-4_Zimmermann_S6_L001_R1_001",
         "p637sDiff6CM_Zimmermann_S46_L005_R1_001",
         "p722s3C190604_Tiburcy_S37_L005_R1_001",
         "p786sC190924A_Tiburcy_S2_L001_R1_001")


#spreaded <- spreaded[, metdat$Sample]

spreaded$Species.name <- rownames(spreaded)
meltd <- melt(spreaded)
head(meltd)

lab <- seq(1:7)

meltd$Species.name <- as.factor(meltd$Species.name)

# !!! Create a custom color scale -- keeps colours constant across plots of the same factors 

myColors <- tableau_color_pal(palette = "Winter")(10)
names(myColors) <- levels(meltd$Species.name)
colScale <- scale_fill_manual(name = "Species.name",values = myColors)


a <- meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity") +
  #scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Number Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  colScale 



b <- meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity", position = "fill") +
  scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Percentage of Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  colScale 

ggsave("filledSpecies.png")

c <- meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~Species.name)+
  theme_minimal()+
  labs(x="Samples",
       y="Number of Reads Mapped")+
  scale_x_discrete(label=lab) +
  guides(fill=FALSE)+
  theme(panel.grid.major = element_blank(), 
        panel.border = element_rect(colour = "grey", 
                                    fill = NA, 
                                    size = 1.0),
        strip.text.x = element_text(size=8),
        axis.text.y = element_text(size = 7)
        ) +
  colScale 
#scale_fill_npg()

ggsave("facetSpecies.png")

library(patchwork)
pat <- ((a + b)/ c) 

pat[[1]] <- pat[[1]] + plot_layout(tag_level = 'new',
      guides = "collect") 
tiff("bacSpecies1.tiff",units="in", width=10,height=7, res=300, compression = 'lzw') 
pat + plot_annotation(tag_levels = c('A','1'),
                      tag_sep = '.',
                      title = 'Possible Bacterial Contaminants at Species Level')

dev.off()


#######################



dat <- na.omit(dat)

head(dat)
dim(dat)
dat$X <- as.factor(dat$X)


new <- melt(dat)


row.names(dat) <- dat$X
dat <- dat[,-1]

head(dat)
samp <- colnames(dat) 

colnames(dat) <- c("X", (seq(1:7)))
colSums(dat[,-1])
dat  <- dat %>%
  filter(X %in% topTen$X)
melted <- melt(dat)

scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)

dev.off()
heatmap(as.matrix(dat[,-1]),Rowv = NA, Colv = NA, col = scaleyellowred)



dat[] <- lapply(dat, function(x) as.numeric(as.character(x)))
pheatmap(as.matrix(dat[,-1]),
         scale = "none",
         show_colnames = FALSE,
         #breaks= max(8),
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="BrBG")))(100))


melted$X <- factor(melted$X,levels = unique(melted$X))
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")
library(ggplot2)
ggplot(melted,aes(x = variable, y = X)) + 
  geom_point(aes(size = value, fill = X), alpha = 0.75, shape = 21)+  
  scale_size_continuous(limits = c(0.000000, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        #legend.title = element_blank(),
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +  
 guides(fill=FALSE)  

############# modified !!!
##############################VIRUS###################
metdat <- read_xlsx("data/basicStats.xlsx")
metdat <- metdat[-4,]
metdat$unmapped <-metdat$`Total Reads` - rowSums(metdat[,-(1:2)])


metdat$Sample <- as.factor(metdat$Sample)
metdat$numbered <- as.factor(metdat$numbered)
f_files<- list.files("data/decon/virus",
                     pattern = "*summary*", full.names = T) 



read_in_feature_counts<- function(file){
  cnt<- read.table(file, sep="\t", header = T)
  cnt <- cnt %>%
    mutate(fileName=file)
  return(cnt)
}


raw_counts <- map(f_files, read_in_feature_counts)

raw_counts_df<- purrr::reduce(raw_counts, union) 

raw_counts_df <- as.data.frame(raw_counts_df)

raw_counts_df$fileName <- gsub("data/decon/virus/unmapped_","\\", raw_counts_df$fileName)
raw_counts_df$fileName <- gsub("_RIBO_UNALIGNED_vs_viruses_subject_summary_all_CT_5.txt","\\", raw_counts_df$fileName)



head(raw_counts_df)
colnames(raw_counts_df)
spreaded <- tidyr::spread(raw_counts_df, fileName, Count)

spreaded[is.na(spreaded)] <- 0


#metdat$Sample %in% colnames(spreaded)
rownames(spreaded) <- spreaded$Species.name
spreaded <- spreaded[,-(1)]
spreaded <- spreaded %>%
  select("SRR1663123_GSM1554465_Fetal_heart_RNA_seq_rep1_Homo_sapiens_RNA-Seq",
         "SRR6706796_GSM2991857_45_HM_3_FC2_Homo_sapiens_RNA-Seq",
         "Sample_r733sCDICM3_hs_totalrna_sr_Zimmermann_D_3",
         "p556sCM10-3-4_Zimmermann_S6_L001_R1_001",
         "p637sDiff6CM_Zimmermann_S46_L005_R1_001",
         "p722s3C190604_Tiburcy_S37_L005_R1_001",
         "p786sC190924A_Tiburcy_S2_L001_R1_001")

## ordering column based on the same order as before 


spreaded$Species.name <- rownames(spreaded)
spreaded <- as.data.frame(spreaded)

meltd <- melt(spreaded)
head(meltd)

lab <- seq(1:7)



# !!! Create a custom color scale -- keeps colours constant across plots of the same factors 
meltd$Species.name <- gsub("NC_[0-9][0-9][0-9][0-9][0-9][0-9].[0-9]","\\", meltd$Species.name)
meltd$Species.name <- as.factor(meltd$Species.name)

myColors <- tableau_color_pal(palette = "Winter")(7)
names(myColors) <- levels(meltd$Species.name)
colScale <- scale_fill_manual(name = "Species.name",values = myColors)
meltd$Species.name <- reorder(meltd$Species.name, meltd$value)
meltd$Species.name <- factor(meltd$Species.name, levels=rev(levels(meltd$Species.name)))


g <- meltd %>%
  mutate(logged=log(value+1, 10)) %>%
  ggplot(aes(x=variable, y=value, fill=Species.name)) +
  geom_bar(stat="identity") +
  scale_x_discrete(label=lab) +
  #scale_y_log10(labels = function(x) format(x,scientific=FALSE))+
  theme_minimal()+
  #scale_y_log10()+
  labs(x="",
       y="Number Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  colScale 
g
##dplyr not in ! 
i <- meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity", position = "fill") +
  scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Percentage of Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  guides(fill=FALSE)+
  colScale 
i

pat <- (g + i)

pat[[1]] <- pat[[1]] + plot_layout(tag_level = 'new',
                                   guides = "collect") 
tiff("00viralSp.tiff",units="in", width=9,height=3, res=300, compression = 'lzw') 
pat + plot_annotation(tag_levels = c('A')) +
  plot_layout(guides="collect")

dev.off()

########################
##############################BACTERIA###################

f<- list.files("data/decon/bacteria",
               pattern = "*summary_ge*", full.names = T)
read_in_feature_counts<- function(file){
  cnt<- read.table(file, sep="\t", header = T)
  cnt <- cnt %>%
    mutate(fileName=file)
  return(cnt)
}

raw_counts <- map(f, read_in_feature_counts)

raw_counts_df<- purrr::reduce(raw_counts, union) 

raw_counts_df <- as.data.frame(raw_counts_df)
raw_counts_df$fileName <- gsub("data/decon/bacteria/unmapped_","\\", raw_counts_df$fileName)
raw_counts_df$fileName <- gsub("_RIBO_UNALIGNED_vs_bacteria_subject_summary_ge_CT_5.txt","\\", raw_counts_df$fileName)

raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum) %>%
  summarize(sum(Count))


raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum) %>%
  arrange(desc(Count)) %>%
  top_n(10) %>%
  summarize(sum(Count))

##87.7 % contributed by first 10 genus
topTen <- raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum) %>%
  arrange(desc(Count)) %>%
  top_n(10) 

metdat$Sample %in% raw_counts_df$fileName

head(raw_counts_df)

raw_counts_df <- raw_counts_df %>%
  filter(Species.name %in% topTen$Species.name)

spreaded <- tidyr::spread(raw_counts_df, fileName, Count)
spreaded[is.na(spreaded)] <- 0

metdat$Sample %in% colnames(spreaded)
rownames(spreaded) <- spreaded$Species.name
spreaded <- spreaded[,-1]
## ordering column based on the same order as before 
spreaded <- spreaded %>%
  select("SRR1663123_GSM1554465_Fetal_heart_RNA_seq_rep1_Homo_sapiens_RNA-Seq",
         "SRR6706796_GSM2991857_45_HM_3_FC2_Homo_sapiens_RNA-Seq",
         "Sample_r733sCDICM3_hs_totalrna_sr_Zimmermann_D_3",
         "p556sCM10-3-4_Zimmermann_S6_L001_R1_001",
         "p637sDiff6CM_Zimmermann_S46_L005_R1_001",
         "p722s3C190604_Tiburcy_S37_L005_R1_001",
         "p786sC190924A_Tiburcy_S2_L001_R1_001")

#spreaded <- spreaded[, metdat$Sample]

spreaded$Species.name <- rownames(spreaded)
meltd <- melt(spreaded)
head(meltd)

lab <- seq(1:7)

meltd$Species.name <- as.factor(meltd$Species.name)

# !!! Create a custom color scale -- keeps colours constant across plots of the same factors 

myColors <- tableau_color_pal(palette = "Winter")(10)
names(myColors) <- levels(meltd$Species.name)
colScale <- scale_fill_manual(name = "Species.name",values = myColors)
meltd$Species.name <- reorder(meltd$Species.name, meltd$value)

x <- meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity") +
  #scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Number Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  colScale
x
y <-  meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity", position = "fill") +
  scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Percentage of Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  colScale 
y


pat <- (x + y)

pat[[1]] <- pat[[1]] + plot_layout(tag_level = 'new',
                                   guides = "collect") 
tiff("00bacGenus1.tiff",units="in", width=9,height=3, res=300, compression = 'lzw') 
pat + plot_annotation(tag_levels = c('A'))+
  plot_layout(guides="collect")

dev.off()


########## species level ##### 

f<- list.files("data/decon/bacteria",
               pattern = "*summary_sp*", full.names = T)
read_in_feature_counts<- function(file){
  cnt<- read.table(file, sep="\t", header = T)
  cnt <- cnt %>%
    mutate(fileName=file)
  return(cnt)
}

raw_counts <- map(f, read_in_feature_counts)

raw_counts_df<- purrr::reduce(raw_counts, union) 


raw_counts_df <- as.data.frame(raw_counts_df)
raw_counts_df$fileName <- gsub("data/decon/bacteria/unmapped_","\\", raw_counts_df$fileName)
raw_counts_df$fileName <- gsub("_RIBO_UNALIGNED_vs_bacteria_subject_summary_sp_CT_5.txt","\\", raw_counts_df$fileName)
dplyr::filter(raw_counts_df, grepl('P', Species.name))

head(raw_counts_df)

raw_counts_df[is.na(raw_counts_df)] <- 0
raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum)  %>%
  summarize(sum(Count))

raw_counts_df %>%
  filter(fileName == "p722s3C190604_Tiburcy_S37_L005_R1_001") %>%
  arrange(desc(Count)) %>%
  select(-fileName)%>%
  top_n(2)%>%
  summarise(sum(Count))


raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum) %>%
  arrange(desc(Count)) %>%
  top_n(20) %>%
  summarize(sum(Count))

##57.7 % contributed by first 10 species
topTen <- raw_counts_df %>%
  group_by(Species.name) %>%
  summarise_at(vars(Count), sum) %>%
  arrange(desc(Count)) %>%
  top_n(10) 

metdat$Sample %in% raw_counts_df$fileName

head(raw_counts_df)

raw_counts_df <- raw_counts_df %>%
  filter(Species.name %in% topTen$Species.name)

spreaded <- tidyr::spread(raw_counts_df, fileName, Count)
spreaded[is.na(spreaded)] <- 0

metdat$Sample %in% colnames(spreaded)
rownames(spreaded) <- spreaded$Species.name

spreaded <- spreaded[,-1]
## ordering column based on the same order as before 
spreaded <- spreaded %>%
  select("SRR1663123_GSM1554465_Fetal_heart_RNA_seq_rep1_Homo_sapiens_RNA-Seq",
         "SRR6706796_GSM2991857_45_HM_3_FC2_Homo_sapiens_RNA-Seq",
         "Sample_r733sCDICM3_hs_totalrna_sr_Zimmermann_D_3",
         "p556sCM10-3-4_Zimmermann_S6_L001_R1_001",
         "p637sDiff6CM_Zimmermann_S46_L005_R1_001",
         "p722s3C190604_Tiburcy_S37_L005_R1_001",
         "p786sC190924A_Tiburcy_S2_L001_R1_001")


#spreaded <- spreaded[, metdat$Sample]

spreaded$Species.name <- rownames(spreaded)
meltd <- melt(spreaded)
head(meltd)

lab <- seq(1:7)

meltd$Species.name <- as.factor(meltd$Species.name)

# !!! Create a custom color scale -- keeps colours constant across plots of the same factors 

myColors <- tableau_color_pal(palette = "Winter")(10)
names(myColors) <- levels(meltd$Species.name)
colScale <- scale_fill_manual(name = "Species.name",values = myColors)
meltd$Species.name <- reorder(meltd$Species.name, meltd$value)


a <- meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity") +
  #scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Number Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  colScale 

a

b <- meltd %>%
  ggplot((aes(x=variable, y=value, fill=Species.name)))+
  geom_bar(stat="identity", position = "fill") +
  scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  scale_x_discrete(label=lab) +
  theme_minimal()+
  labs(x="",
       y="Percentage of Reads Mapped")+
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  colScale 
b


pat <- (a + b)

pat[[1]] <- pat[[1]] + plot_layout(tag_level = 'new',
                                   guides = "collect") 
tiff("00bacSpecies1.tiff",units="in", width=9,height=3, res=300, compression = 'lzw') 
pat + plot_annotation(tag_levels = c('A'))+
  plot_layout(guides="collect")

dev.off()

pat <- (g + h + i) / (x + y)/( a + b)
pat[[1]] <- pat[[1]] + plot_layout(tag_level = 'new',
                                   guides = "collect")
pat[[2]] <- pat[[2]] + plot_layout(tag_level = 'new',
                                   guides = "collect")
pat[[3]] <- pat[[3]] + plot_layout(tag_level = 'new',
                                   guides = "collect")
pat + plot_annotation(tag_levels = c('A','1'),
                      tag_sep = '.')



