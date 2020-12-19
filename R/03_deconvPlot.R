## all the deconvolution plots from fig 4.11 to 4.15


cib <- read.csv("/Users/bhuvaneswari/Downloads/CIBERSORTx_Job30_Results.csv")

ci1 <- cib %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "groups", values_to = "proportions")

ci2 <- dcast(ci1, groups~type, mean, na.rm=TRUE, value.var = "Correlation")
ci3 <- ci2 %>%
  pivot_longer(c('Adult','Atrial CM','CM','Fetal','Ventricular CM','EHM'), 
               names_to = "categ",
               values_to = "Correlation")

ci4 <- ci3[1:6,]
ci5 <- ci3[1:6,]
mean(ci5$RMSE)
mean(ci4$corr)
###########
cib1 <- cib %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "groups", values_to = "proportions")

cib2 <- dcast(cib1, type~groups, mean, na.rm=TRUE, value.var = "proportions")
cib3 <- cib2 %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "well", values_to = "proportions")
cib4 <- cib3 %>%
  group_by(type) %>%
  mutate(percent=(round(100*proportions/sum(proportions),0)))

#write_csv(cib6,"deconProp1.csv")
cib6 <- cib4 %>%
  filter(str_detect(Mixture,'GMP'))
cib4$Mixture


library(readr)
modCib <- read.csv("data/deconProp.csv")
tiff("03_02Deconv-MDC.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
#create the stacked bar plot based on your data
modCib$well <- factor(modCib$well, levels=c("Comparability Factor","d15.S1","d30.S1","d15.S2","d30.S2"))
barPlot <- ggplot(data = modCib, aes(y=proportions, x=type, fill=well)) + 
  geom_bar(stat="identity", width = 0.5, position="fill") + 
  xlab('') + ylab('') +
  #use JOELS great solution for the label position 
  #and add percentage based on variable 'relative', otherwise use 'number'
  geom_text(aes(x = type, label = paste0(percent,'%')),
            colour = 'grey10',position=position_fill(vjust=0.5), size=6) + 
  labs(fill='') + theme_bw() +
  scale_fill_manual(values=c("d15.S1"="#A0CBE8",
                             "d15.S2" ="#FFBE7D",
                             "d30.S1"="#4E79A7",
                             "d30.S2"="#F28E2B",
                             "Comparability Factor" = "grey"))+
  #scale_fill_tableau(palette="Tableau 20")+
  theme_few()+
  theme(
    axis.text.y.left = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 15, face="bold")
  )
barPlot
dev.off()
##############
cib1 <- cib %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "groups", values_to = "proportions")

cib2 <- dcast(cib1, Mixture~groups, mean, na.rm=TRUE, value.var = "proportions")
cib3 <- cib2 %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "well", values_to = "proportions")
cib4 <- cib3 %>%
  group_by(Mixture) %>%
  mutate(percent=(round(100*proportions/sum(proportions),0)))
library(stringr)
cib6 <- cib4 %>%
  filter(str_detect(Mixture,'GMP'))
cib4$Mixture

cib7 <- cib4 %>%
  filter(str_detect(Mixture,'Adult_Heart'))
library(numform)
cib6$mixture <- as.factor(cib6$mixture)
boxCM <- ggplot(cib6, aes( well, percent)) +
  geom_boxplot(aes(fill = stage(well, after_scale = alpha(fill, 0.9))))+
  scale_fill_manual(values=c("d15.S1"="#A0CBE8",
                              "d15.S2" ="#FFBE7D",
                              "d30.S1"="#4E79A7",
                              "d30.S2"="#F28E2B"))+
  theme_few()+
  guides(fill=FALSE)+
  scale_y_continuous(labels= ff_percent(digits=0))+
  theme(
    axis.title=element_blank(),
    axis.text.y=element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, face="bold", angle=45, hjust=1)
  )
#####################################################
tiff("02_02Deconv-MDC.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
#create the stacked bar plot based on your data

cib6$well <- factor(cib6$well, levels=c("d15.S1","d30.S1","d15.S2","d30.S2"))
barCM <- ggplot(data = cib6, aes(y=proportions, x=Mixture, fill=well)) + 
  geom_bar(stat="identity", width = 0.5, position="fill") + 
  xlab('') + ylab('') +
  #use JOELS great solution for the label position 
  #and add percentage based on variable 'relative', otherwise use 'number'
  geom_text(aes(x = Mixture, label = paste0(percent,'%')),
            colour = 'grey10',position=position_fill(vjust=0.5), size=3) + 
  labs(fill='') + theme_bw() +
  scale_fill_manual(values=c("d15.S1"="#A0CBE8",
                             "d15.S2" ="#FFBE7D",
                             "d30.S1"="#4E79A7",
                             "d30.S2"="#F28E2B"))+
  #scale_fill_tableau(palette="Tableau 20")+
  theme_few()+
  theme(
    axis.text.y.left = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, face="bold", angle=45, hjust=1)
  )

tiff("01deconvCM.tiff",units="in", width=6,height=3, res=300, compression = 'lzw') 
boxCM + barCM +
  plot_annotation(tag_levels = "A")
dev.off()
###############

library(gt)

deconvTab <- tribble(
  ~"Day 15", ~"Day 30",
  "Non-Contractile (d15:S1)", "Non-Contractile (d30:S1)",
  "Cardiac-like Cells (d15:S2)","Cardiomyocytes (d30:S2)"
)
gt_tbl <- gt(data=deconvTab)

table <- gt_tbl %>%
  tab_header(
    title = "2 Time Points & their Cell Populations",
    subtitle = "As per single cell reference data"
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#A0CBE8")
    ),
    locations = cells_body(
      columns = vars("Day 15"),
      rows = 1)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#FFBE7D")
    ),
    locations = cells_body(
      columns = vars("Day 15"),
      rows = 2)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#4E79A7")
    ),
    locations = cells_body(
      columns = vars("Day 30"),
      row = 1)
  )%>%
  tab_style(
    style = list(
      cell_fill(color = "#F28E2B")
    ),
    locations = cells_body(
      columns = vars("Day 30"),
      row = 2)
  )%>%
  cols_align(
    align = "center"
  )
library(png)
library(patchwork)
library(ggplot2)
library(gtable)
library(magick)

img1 <- image_read('data/graphAbs.png')
tiff("deconvBar.tiff",units="in", width=6,height=5, res=300, compression = 'lzw') 
barPlot
dev.off()
tablegrid.arrange(barPlot, image, byrow = TRUE) 
###########################################
c <- read.csv("data/CIBERSORTx_Job37_Results.csv")
c
c1 <- c %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "groups", values_to = "proportions")

c2 <- dcast(c1, Mixture~groups, mean, na.rm=TRUE, value.var = "proportions")
c3 <- c2 %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "well", values_to = "proportions")
c4 <- c3 %>%
  group_by(Mixture) %>%
  mutate(percent=(round(100*proportions/sum(proportions),0)))

write_csv(cib6,"deconProp1.csv")
c6 <- c4 %>%
  filter(str_detect(Mixture,'EHM'))
c6$Mixture
c6 <- c6 %>%
  mutate(fromWhere = case_when(Mixture = (str_detect(Mixture,'medium')) ~ 'PRJNA362579',
                               TRUE ~ 'In-House'),
         colorCode = case_when(fromWhere == "PRJNA362579" ~ "black",
                               TRUE ~"red"))
colours <- c6$colorCode
c6$fromWhere <- as.factor(c6$fromWhere)
library(readr)
modCib <- read.csv("deconProp.csv")
tiff("02_02Deconv-MDC.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
#create the stacked bar plot based on your data
c6$well <- factor(c6$well, levels=c("d15.S1","d30.S1","d15.S2","d30.S2"))
filtout <- c6 %>%
  filter(fromWhere == "In-House")
filtout$well
down1 <- ggplot(data = filtout, aes(y=proportions, x=Mixture, fill=well)) + 
  #facet_wrap(~Mixture)+
  geom_bar(stat="identity", width = 0.5, position="fill") + 
  xlab('') + ylab('') +
  #use JOELS great solution for the label position 
  #and add percentage based on variable 'relative', otherwise use 'number'
  geom_text(aes(x = Mixture, label = paste0(percent,'%')),
            colour = 'grey10',position=position_fill(vjust=0.5), size=3) + 
  labs(fill='') + theme_bw() +
  scale_fill_manual(values=c("d15.S1"="#A0CBE8",
                             "d15.S2" ="#FFBE7D",
                             "d30.S1"="#4E79A7",
                             "d30.S2"="#F28E2B",
                             "Comparability Factor" = "grey"))+
  #scale_fill_tableau(palette="Tableau 20")+
  theme_few()+
  theme(
    axis.text.y.left = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.x = element_text(size = 10,hjust=1, face="bold", angle=45, 
    #                           color =c(rep("#cf3e53",10),rep("#6b6b6b",7))), 
    axis.text.x = element_text(size=10, hjust=1, angle=45)
  )
down1
okay <- element_text(color=c6$colorCode)
up <- ggplot(c6, aes( well, percent)) +
  geom_boxplot(aes(fill = stage(well, after_scale = alpha(fill, 0.9))))+
  facet_wrap(~fromWhere)+
  scale_fill_manual(values=c("d15.S1"="#A0CBE8",
                             "d15.S2" ="#FFBE7D",
                             "d30.S1"="#4E79A7",
                             "d30.S2"="#F28E2B"))+
  theme_few()+
  guides(fill=FALSE)+
  scale_y_continuous(labels= ff_percent(digits=0))+
  theme(
    axis.title=element_blank(),
    axis.text.y=element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, face="bold", angle=45, hjust=1)
  )
dev.off()

tiff("01deconvEHM.tiff",units="in", width=13,height=10, res=300, compression = 'lzw') 

down <- down1 | down2
pa <- up / down
pa[[2]] <- pa[[2]] + plot_layout(tag_level = 'new',
                                 guides = "collect") 
pa + 
  plot_annotation(tag_levels = c('A','1'),
                  tag_sep = '.')
#plot_layout(guides = "collect")
dev.off()
##############################
corr <- c %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "groups", values_to = "proportions")

corr <- dcast(corr, Mixture~groups, mean, na.rm=TRUE, value.var = "Correlation")
corr <- corr %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "well", values_to = "Correlation")
corr <- corr %>%
  mutate(whichGroup = case_when(Mixture = (str_detect(Mixture,'Adult_Heart')) ~ 'Adult Heart',
                                Mixture = (str_detect(Mixture,'CM')) ~ 'CM',
                                Mixture = (str_detect(Mixture,'EHM')) ~ 'EHM',
                                Mixture = (str_detect(Mixture,'Fetal_Heart')) ~ 'Fetal Heart',
                                Mixture = (str_detect(Mixture,'Fib')) ~ 'Fibroblasts',
                                Mixture = (str_detect(Mixture,'ipsc')) ~ 'iPSC',
                                TRUE ~ 'In-House'))

corred <- corr %>%
  dplyr::select(-Mixture) %>% 
  group_by(whichGroup, well) %>%
  summarise_each(funs(median)) %>%
  arrange(desc(Correlation))

rmse <- c %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "groups", values_to = "proportions")

rmse <- dcast(rmse, Mixture~groups, mean, na.rm=TRUE, value.var = "RMSE")
rmse <- rmse %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "well", values_to = "RMSE")
rmse <- rmse %>%
  mutate(whichGroup = case_when(Mixture = (str_detect(Mixture,'Adult_Heart')) ~ 'Adult Heart',
                                Mixture = (str_detect(Mixture,'CM')) ~ 'CM',
                                Mixture = (str_detect(Mixture,'EHM')) ~ 'EHM',
                                Mixture = (str_detect(Mixture,'Fetal_Heart')) ~ 'Fetal Heart',
                                Mixture = (str_detect(Mixture,'Fib')) ~ 'Fibroblasts',
                                Mixture = (str_detect(Mixture,'ipsc')) ~ 'iPSC',
                                TRUE ~ 'In-House'))



rmsed <- rmse %>%
  dplyr::select(-Mixture) %>% 
  group_by(whichGroup, well) %>%
  filter(RMSE== min(RMSE)) %>% 
  arrange(RMSE)


c8 <- c4 %>%
  mutate(whichGroup = case_when(Mixture = (str_detect(Mixture,'Adult_Heart')) ~ 'Adult Heart',
                                Mixture = (str_detect(Mixture,'CM')) ~ 'CM',
                                Mixture = (str_detect(Mixture,'EHM')) ~ 'EHM',
                                Mixture = (str_detect(Mixture,'Fetal_Heart')) ~ 'Fetal Heart',
                                Mixture = (str_detect(Mixture,'Fib')) ~ 'Fibroblasts',
                                Mixture = (str_detect(Mixture,'ipsc')) ~ 'iPSC',
                               TRUE ~ 'In-House'))

c8$whichGroup <- factor(c8$whichGroup)
c8$well <- factor(c8$well)
proped <- c8 %>%
  dplyr::select(-Mixture) %>% 
  group_by(whichGroup, well) %>%
  dplyr::select(-Mixture) %>% 
  summarise_each(funs(median))

grpd <- merge(proped, corred)
grpd <- merge(grpd, rmsed)
grpd$well <- factor(grpd$well, levels=c("d15.S1","d30.S1","d15.S2","d30.S2"))
grpd$whichGroup <- factor(grpd$whichGroup, levels=c("iPSC","Fibroblasts","Adult Heart","CM","EHM","Fetal Heart"))
tiff("groupdDecon2.tiff",units="in", width=6,height=5, res=300, compression = 'lzw') 
one <- ggplot(data = grpd, aes(y=proportions, x=whichGroup, fill=well)) + 
  #facet_wrap(~Mixture)+
  geom_bar(stat="identity", width = 0.5, position="fill") + 
  xlab('') + ylab('') +
  #use JOELS great solution for the label position 
  #and add percentage based on variable 'relative', otherwise use 'number'
  geom_text(aes(x = whichGroup, label = paste0(percent,'%')),
            colour = 'grey10',position=position_fill(vjust=0.5), size=3) + 
  labs(fill='') + theme_bw() +
  scale_fill_manual(values=c("d15.S1"="#A0CBE8",
                             "d15.S2" ="#FFBE7D",
                             "d30.S1"="#4E79A7",
                             "d30.S2"="#F28E2B"))+
  #scale_fill_tableau(palette="Tableau 20")+
  theme_few()+
  theme(
    axis.text.y.left = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.x = element_text(size = 10,hjust=1, face="bold", angle=45, 
    #                           color =c(rep("#cf3e53",10),rep("#6b6b6b",7))), 
    axis.text.x = element_text(size=10, hjust=1, angle=45)
  )

two <- grpd %>% 
  select(whichGroup, RMSE) %>% 
  group_by(whichGroup) %>%
  summarise_each(funs(median)) %>% 
  ggplot() +
  geom_bar(aes(x= whichGroup, y= RMSE, fill=RMSE),width = 0.5, stat="identity") +
  theme_few()+
  theme(
    axis.text.x = element_text(size=10, hjust=1, angle=45)
  ) +
  labs(x=" ", y= " ")+
  scale_fill_continuous_tableau(palette = "Orange")


one / two +
  plot_layout(heights  = c(2, 1)) +
  plot_annotation(tag_levels = "A")
dev.off()
tiff("coord.tiff",units="in", width=10,height=2, res=300, compression = 'lzw') 
ggplot(grpd, aes(x=Correlation, y=c(0), colour=Correlation)) +
  geom_line(size=2) +
  geom_point(size=4) +
  #scale_color_gradient(low="darkblue",high="darkred")+
  scale_colour_viridis_c()+
  geom_text(aes(x = Correlation, label = whichGroup),
            colour = 'grey10',hjust=1,vjust=2, size=4,
            check_overlap = T,
            angle=45) +  
  geom_text(aes(x = Correlation,y=c(0), 
                label = (round(Correlation,2))),
                                   colour = 'grey10',
            hjust=0,
            vjust=-1, 
            size=4,
            check_overlap = T,
            angle=0)+
  theme_void()
dev.off()
tiff("rmsed.tiff",units="in", width=11,height=2, res=300, compression = 'lzw') 
ggplot(grpd, aes(x=RMSE, y=c(0), colour=RMSE)) +
  geom_line(size=2) +
  geom_point(size=4) +
  #scale_color_gradient(low="darkblue",high="darkred")
  geom_text(aes(x = RMSE, label = whichGroup),
            colour = 'grey10',hjust=1,vjust=2, size=4,
            check_overlap = T,
            angle=45) +
  geom_text(aes(x = RMSE,y=c(0), 
                label = (round(RMSE,2))),
            colour = 'grey10',
            hjust=0,
            vjust=-1, 
            size=4,
            check_overlap = T,
            angle=0)+
  scale_x_reverse()+
  scale_colour_viridis_c(direction = -1)+
  theme_void()
dev.off()


  


hmm <- c("red")
##############3
sigMat <- read.table("data/CIBERSORTx_Job24_combined_d15-d30_inferred_phenoclasses.CIBERSORTx_Job24_combined_d15-d30_inferred_refsample.bm.K999.txt",sep="\t", header  = T)
d30.sig2 <- sigMat %>%
  select(NAME, d30.S2) %>%
  
  arrange(desc(d30.S2)) %>%
  slice(1:100)

d15.sig2 <- sigMat %>%
  select(NAME, d15.S2) %>%
  arrange(desc(d15.S2)) %>%
  slice(1:100)

d15.sig1 <- sigMat %>%
  select(NAME, d15.S1) %>%
  arrange(desc(d15.S1)) %>%
  slice(1:100)

d30.sig1 <- sigMat %>%
  select(NAME, d30.S1) %>%
  arrange(desc(d30.S1)) %>%
  slice(1:100)
sigMatGenes <- cbind2(d30.sig2, d30.sig1, d15.sig2, d15.sig2)
sigMatGenes <- cbind(sigMatGenes, d15.sig1)
sigMatGenes <- cbind(sigMatGenes, d15.sig2)
write.csv(sigMatGenes, "sorted100GenesSigMatScRNA.csv")
addToDF <- function(df, v){
  nRow <- nrow(df)
  lngth <- length(v)
  if(nRow > lngth){
    length(v) <- nRow
  }else if(nRow < lngth){
    df[(nRow+1):lngth, ] <- NA
  }
  cbind(df,v)
}

addToDF(sortedGenesSigMat, d15.sig1$NAME)

sortedGenesSigMat <- cbind2(d30.sig2, d30.sig1, d15.sig2, d15.sig2)
sortedGenesSigMat <- addToDF(sortedGenesSigMat, d15.sig1$NAME)
sortedGenesSigMat <- addToDF(sortedGenesSigMat, d15.sig1$d15.S1)
sortedGenesSigMat <- addToDF(sortedGenesSigMat, d15.sig2$NAME)
sortedGenesSigMat <- addToDF(sortedGenesSigMat, d15.sig2$d15.S2)

write.csv(sortedGenesSigMat, "sortedGenesSigMatScRNA.csv")

# Pairwise correlation between samples (columns)
cols.cor <- cor(filtHeat, 
                use = "pairwise.complete.obs", 
                method = "pearson")

# Pairwise correlation between rows (genes)
rows.cor <- cor(t(filtHeat), 
                use = "pairwise.complete.obs", 
                method = "pearson")




library(pheatmap)
mew <- mew[apply(mew, MARGIN = 1, FUN = function(x) sd(x) != 0),]
tiff("heatMapThesis.tiff",units="in", width=10,height=10, res=300, compression = 'lzw') 
pheatmap(
  filtHeat,
  #scale = "row",
  #cluster_cols=hc.complete,
  #cluster_rows =hc.complete.rows,
  clustering_distance_cols = as.dist(1 - cols.cor),
  clustering_distance_rows = as.dist(1 - rows.cor),
  #cutree_cols = 4,
  #cutree_rows = 6,
  annotation_colors = heatmapColScale,
  #annotation_col = colAn,
  #annotation_row = rowAn,
  show_rownames = FALSE,
  #show_colnames = FALSE,
  labels_col = hc.complete$label,
  treeheight_row = 0,
  angle_col = 315,
  fontsize_col = 7
  #col=inferno(100)
  #  col = inferno(length(breaksList)),
  # breaks = breaksList
)

dev.off()

filtered <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  #dplyr::filter(gene %in% top_genes) %>%
  select(gene, PC1, PC4) %>%
  dplyr::filter(PC4 > -0.01 & PC4 < 0.01 ) 

onlyEHM <- dds[,dds$group == 'EHM']
onlyEHM <- vst(onlyEHM)
filtHeat <- assay(onlyEHM) %>%
  as_tibble(rownames = "gene") %>% 
  dplyr::filter(gene %in% sigMat$NAME)

dim(filtHeat)
filtHeat <- as.data.frame(filtHeat)
rownames(filtHeat) <- filtHeat$gene
filtHeat <- filtHeat[,-1]

filtHeat <- filtHeat[apply(filtHeat, MARGIN = 1, FUN = function(x) sd(x) != 0),]

pheatmap(na.omit(filtHeat),
         scale="row",
         method = "complete",
         color = inferno(100),
         cutree_rows = 4,
         cutree_cols = 3,
         show_rownames = FALSE,
         labels_col = paste(onlyEHM$group, onlyEHM$specification, sep = " - " ))

onlyEHM@colData
write.table(as.data.frame(assay(onlyEHM)), "onlyEHMrawCounts.txt")
write.table(as.data.frame(onlyEHM@colData), "onlyEHMmetaData.txt")
onlyEHMmet<- read.table("onlyEHMmetaData.txt")
onlyEHMcnt <- read.table("onlyEHMrawCounts.txt")

dds <- DESeqDataSetFromMatrix(countData = onlyEHMcnt, colData = onlyEHMmet, design=~ fromWhere)
dds <- dds[ rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds)

contrast_oe <- c("fromWhere", "PRJNA362579", "Malte")

res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)
results(dds)
res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken)
res_tableOE %>% data.frame() %>% View()
summary(res_tableOE)
padj.cutoff <- 0.05
lfc.cutoff <- 1
library(tibble)
res_tableOE_tb <- res_tableOE %>%
  as_tibble(rownames ="gene") 
# rownames_to_column(var="gene"
sigOE <- res_tableOE_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)


interesected <- intersect(sigOE$gene, sigMat$NAME)

write(interesected, "checkTheseOut.txt")

############### Deconvolution Pseudo Bulk 

cib <- read.csv("/Users/bhuvaneswari/Downloads/CIBERSORTx_Job30_Results.csv")

ci1 <- cib %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "groups", values_to = "proportions")

ci2 <- dcast(ci1, groups~type, mean, na.rm=TRUE, value.var = "Correlation")
ci3 <- ci2 %>%
  pivot_longer(c('Adult','Atrial CM','CM','Fetal','Ventricular CM','EHM'), 
               names_to = "categ",
               values_to = "Correlation")

ci4 <- ci3[1:6,]
ci5 <- ci3[1:6,]
mean(ci5$RMSE)
mean(ci4$corr)
###########
library(readr)
psdB <- read.table("data/CIBERSORTx_Job47_Results.txt", sep = "\t", header = T)
psdB$renamed <- c("1-P", "2-P", "3-P")
psdB <- psdB %>%
  pivot_longer(c('d30.S2','d30.S1','d15.S1','d15.S2'), names_to = "groups", values_to = "proportions")

psdB <- psdB %>%
  group_by(Mixture) %>%
  mutate(percent=(round(100*proportions/sum(proportions),0)))
write.csv(psdB,"pseudoBA.csv")

tiff("Deconv-pseudoBulk.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
curPM <- read.csv("pseudoBAMod.csv")

ggplot()+
  geom_bar(data = curPM,
             aes(y=percent, x=Renamed, fill=groups),  
           stat="identity", position = "dodge")+
  theme_few() +
  labs(y="Proportions of Sub-Groups (%)", x= "Pseudo Bulk Samples ")+
  scale_fill_manual(values=c("d15.S1"="palegreen",
                              "d15.S2" ="paleturquoise1",
                              "d30.S1"="palevioletred1",
                              "d30.S2"="peachpuff",
                              "d15.S1-P"= "palegreen3",
                              "d15.S2-P" ="paleturquoise3",
                              "d30.S1-P"="palevioletred3",
                              "d30.S2-P"="peachpuff3")) +
  theme(legend.title = element_blank())
dev.off()
