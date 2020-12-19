library(readxl)
library(readr)
library(ggplot2)
library(reshape2)
library(viridis)
library(ggthemes)
library(ggforce)
library(purrr)
library(tidyverse)
library(DESeq2)
library(readxl)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ggsci)
library(PCAtools)
library(EBSeq)
library(patchwork)
library(pheatmap)
library(rafalib)
library(DESeq2)
library(pvclust)
library(rafalib)
library(PCAtools)
library(purrr)
library(readr)
library(dplyr)
library(readxl)
library(ggbeeswarm)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(wesanderson)
library(gridExtra)
library(grid)
library(gghighlight)
library(plotly)
library(broom)
library(corrplot)
library(pheatmap)
library(factoextra)
############
heartGenes <- read_tsv("data/tissue_category_rna_heart-2.tsv")
heartGenes <- heartGenes$Gene
heartGenesOne <- read_tsv("data/tissue_category_rna_heart.tsv")
heartGenesTwo <- read_tsv("data/tissue_category_rna_heart-3.tsv")
heartGenesEnriched <- heartGenesOne$Gene
heartGenesFound <- heartGenesTwo$Gene

gen <- unique(c(heartGenes,heartGenesFound,heartGenesFound))

############
cnt <- read.table("data/rawCountsThesis.txt")

metdata <- read_csv("data/metaDataThesis.csv")
colnames(cnt) == metdata$x
metdata <- as.data.frame(metdata)
row.names(metdata) <- metdata$x
colnames(metdata)
metdata %>%
  filter(str_detect(group,'EHM'))

dds <- DESeqDataSetFromMatrix(countData = cnt, colData = metdata, design=~ Project_AccessionNumber)
dds <- dds[ rowSums(counts(dds)) > 5, ]
forDeconv <- assay(dds)
colnames(forDeconv) <- paste(dds$group, dds$specification, sep = " - " )
forDeconv
write.csv(as.data.frame(forDeconv), "forDeconv.csv", row.names =TRUE,quote=FALSE)
?write_csv
vsd <- vst(dds)
vsd$group <- as.factor(vsd$group)

vsd <- vsd[,vsd$group != "ipsc"]
vsd <- vsd[,vsd$group != "Fib"]
vsd <- vsd[,vsd$group != "Rh"]
vsd <- vsd[,vsd$fromWhere != "Cyg"]


colscale <- c("Adult_Heart" = "#e15759",
  "Fetal_Heart" = "#9d983d",
  "EHM"="#90728f",
  "CM"="#f3a546")


vsd$group <- fct_relevel(vsd$group,"Adult_Heart",
                         "Fetal_Heart",
                         "EHM",
                         "CM")

# perform pca using prcomp func, choosing the top n number of genes 
ntop=2000
Pvars <- rowVars(assay(vsd))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]

PCA <- prcomp(t(assay(vsd)[select, ]), scale = T)
okay <- vsd[,vsd$group == "EHM"]

x <- as.data.frame(l)
x$samp <- rownames(x)
x <- x %>% 
  filter(samp %in% okay$x)
x <- x %>% 
  select(1:4, samp)
x$try <- okay$fromWhere
library(reshape2)
library(ggplot2)
library(ggridges)

d <- melt(x, id = 1:4, measure = 1:ncol(x))
d <- melt(x)

d <- d %>%
  mutate(source = case_when(try == "Malte" ~ 'In-House',
                               TRUE ~ 'PRJNA362579'))

ggplot(as.data.frame(d), aes(x=value, fill=source)) + 
  geom_density(alpha=0.5) +
  facet_wrap(~variable) +
  labs(y=" ")+
  scale_fill_wsj() +
  theme_linedraw()+
  theme(panel.background = element_blank()) +
  theme(strip.background = element_rect(fill="beige"))+
  theme(strip.text = element_text(colour = 'black')) +
  guides(fill=FALSE)

tiff("facetWrapEHMs",units="in", width=5,height=4, res=300, compression = 'lzw') 
ggplot(as.data.frame(d), aes(x=value, y=source, fill=source)) + 
  geom_density_ridges2(alpha=0.5, scale=2,
                       jittered_points = TRUE,
                       position = position_points_jitter(width = 0.05, height = 0),
                       point_shape = '|', point_size = 3, point_alpha = 1) +
  facet_wrap(~variable)+
 # theme_few()+
  labs(y=" ",
       x=" ")+
  scale_fill_wsj() +
  theme_linedraw()+
  theme(panel.background = element_blank()) +
  theme(strip.background = element_rect(fill="beige"))+
  theme(strip.text = element_text(colour = 'black')) +
  guides(fill=FALSE)
dev.off()



ggplot(as.data.frame(l), aes(x=PC1, fill=vsd$group)) + 
  geom_density(alpha=0.5)
  

l <- PCA$x
#l <- plotPCA(vsd, intgroup=c("group","specification"), returnData=TRUE)
#percent variation each component accounts for
percentVarp <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

ggplot(as.data.frame(l), aes(PC2, PC4, color = vsd$group))+
  geom_point(size =3) +
  facet_wrap(~vsd$group)+
  #ggtitle("PCA-using all genes") +
  theme(plot.title = element_text(hjust=0, size=8))+
  #geom_label_repel(data = as.data.frame(l),aes(PC1, PC4),
  #             label = vsd$specification) +
  xlab(paste0("PC2: ", round(percentVarp[2],digits=2),"% variance"))+
  ylab(paste0("PC4: ", round(percentVarp[4],digits = 2), "% variance"))+
  guides(shape=FALSE)+
  scale_color_manual(values = colscale)+
  theme_few()+
  theme(legend.title = element_blank())

one <- ggplot(as.data.frame(l), aes(PC1, PC2, color = vsd$group))+
  geom_point(size =3) +
  #ggtitle("PCA-using all genes") +
  theme(plot.title = element_text(hjust=0, size=8))+
  #geom_label_repel(data = as.data.frame(l),aes(PC1, PC4),
   #             label = vsd$specification) +
  xlab(paste0("PC1: ", round(percentVarp[1],digits=2),"% variance"))+
  ylab(paste0("PC2: ", round(percentVarp[2],digits = 2), "% variance"))+
  guides(shape=FALSE)+
 scale_color_manual(values = colscale)+
  theme_few()+
  theme(legend.title = element_blank())

two <- ggplot(as.data.frame(l), aes(PC1, PC4, color = vsd$group))+
  geom_point(size =3) +
  #ggtitle("PCA-using all genes") +
  theme(plot.title = element_text(hjust=0, size=8))+
  #geom_label_repel(data = as.data.frame(l),aes(PC1, PC4),
  #             label = vsd$specification) +
  xlab(paste0("PC1: ", round(percentVarp[1],digits=2),"% variance"))+
  ylab(paste0("PC4: ", round(percentVarp[4],digits = 2), "% variance"))+
  guides(shape=FALSE)+
  scale_color_manual(values = colscale)+
  theme_few()+
  theme(legend.title = element_blank())

tiff("twoPCA.tiff",units="in", width=9, height=3,res=300, compression = 'lzw') 
(one | two) + plot_layout(guides="collect") + plot_annotation(tag_levels = "A")
dev.off()
ordOne <- ggplot(as.data.frame(l), aes(vsd$group, PC1, fill = vsd$group))+
  geom_boxplot()+
  scale_fill_manual(values=colscale) +
  xlab(" ")+
  theme_few()+
  theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "bold"),
        legend.title = element_blank())

ordTwo <- ggplot(as.data.frame(l), aes(vsd$group, PC2, fill = vsd$group))+
  geom_boxplot() +
  scale_fill_manual(values=colscale) +
  xlab(" ")+
  theme_few()+
  theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "bold"),
        legend.title = element_blank())

  
library(ggplotify)
screePlot <- ggplot(df.eig, aes(dimension, cumulative.variance.percent))+
  geom_bar(stat = "identity", fill="beige", color="grey")+ 
  geom_line(color = "grey") + 
  geom_point(shape = 19, color = "red")+
  geom_text(label = text_labels$cumulative.variance.percent, vjust = -0.8, 
            hjust = 0.8)+
  scale_x_continuous(breaks = c(1:8))+
  scale_y_continuous(limits = c(0, 100), breaks=seq(0.00,100.00,by=20))+
  labs(y="Percentage of Cummulative Variance Explained",
       x="Dimension")+
  theme_few()+
  theme(axis.text.x = element_text(face="bold"))

tiff("modPCAtop.tiff",units="in", width=12, height=4.5,res=300, compression = 'lzw') 
(screePlot | ordOne | ordTwo ) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c("A","1"))
dev.off()
## trying a combination of PCA and loadings plot


######## blah

two <- ggplot(as.data.frame(l), aes(PC1, PC4, color = vsd$group))+
  geom_point(size =3) +
  #ggtitle("PCA-using all genes") +
  theme(plot.title = element_text(hjust=0, size=8))+
  #geom_label_repel(data = as.data.frame(l),aes(PC1, PC4),
  #             label = vsd$specification) +
  xlab(paste0("PC1: ", round(percentVarp[1],digits=2),"% variance"))+
  ylab(paste0("PC4: ", round(percentVarp[4],digits = 2), "% variance"))+
  guides(shape=FALSE)+
  scale_color_manual(values = colscale)+
  theme_few()+
  theme(legend.title = element_blank())

# extract values
PCA$rotation %>% 
  # convert to a tibble
  as_tibble(rownames = "gene")


tidy(PCA, matrix = "variables")
top_genes <- PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "1" |  PC == "2") %>%
  # for each PC
  group_by(PC) %>%
  # sort descending value
  dplyr::arrange(desc(abs(value))) %>%
  # take top 5 rows of each PC group
  dplyr::slice(1:100) %>%
  # extract the column (gene name) from the table
  pull(column) %>%
  # retain unique gene names only
  unique()

ok <- PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "4") %>%
  dplyr::arrange(desc(abs(value))) %>% 
  mutate(mod = value^2*100) 
  
sum(ok$mod)

top_genes <- PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "1") %>%
  # sort descending value
  dplyr::arrange(desc(abs(value))) 

which(top_genes$column=="GABRA4")

write.csv(top_genes,"geneLoadingPC3.csv")

gene_loadings <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::filter(gene %in% top_genes) 

PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::select(gene, PC4) %>% 
  dplyr::arrange(desc(abs(PC4))) %>%
  mutate(mut=(PC4^2*100))

####### Making the table
p1 <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::select(gene, PC1) %>% 
  dplyr::arrange(desc(abs(PC1))) %>%
  slice(1:50)

p2 <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::select(gene, PC2) %>% 
  dplyr::arrange(desc(abs(PC2))) %>%
  slice(1:50)

p3 <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::select(gene, PC3) %>% 
  dplyr::arrange(desc(abs(PC3))) %>%
  slice(1:50)

p4 <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::select(gene, PC4) %>% 
  dplyr::arrange(desc(abs(PC4))) %>%
  slice(1:50)

bound <- cbind(p1,p2,p3,p4)
write.csv(bound, "results/pcAndGenes.csv")
######

filtered <- gene_loadings %>%
  dplyr::filter(gene %in% gen)%>%
  select(gene, PC1, PC2)


filtered <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::filter(gene %in% heartGenes) %>%
  select(gene, PC1, PC2)
  dplyr::filter(PC2 > -0.01 & PC2 < 0.01 )  

##


###

ggplot() +
  geom_segment(data=gene_loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),alpha=0.1,  
               arrow = arrow(length = unit(0.1, "in")),
               colour = "#e15759",
               linetype=1) +
  geom_text_repel(data = filtered, 
                   aes(x = PC1, y = PC2, label = gene),
                  size = 3,
                  segment.size = 0.5,
                  segment.alpha = 0.4,
                  segment.colour = "#cccccc",
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_few()


PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "1" |  PC == "2")


##### loadings plot for PC1 vs PC4 
top_genes4 <- PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "4"| PC == "1") %>%
  # for each PC
  group_by(PC) %>%
  # sort descending value
  dplyr::arrange(desc(abs(value))) %>%
  # take top 5 rows of each PC group
  dplyr::slice(1:50) %>%
  # extract the column (gene name) from the table
  pull(column) %>%
  # retain unique gene names only
  unique()

gene_loadings4 <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::filter(gene %in% top_genes4) 


ggplot(gene_loadings4) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC4),alpha=0.1,  
               arrow = arrow(length = unit(0.1, "in")),
               colour = "#e15759",
               linetype=1) +
  geom_text_repel(data = gene_loadings4, 
                  aes(x = PC1, y = PC4, label = gene),
                  size = 3,
                  segment.size = 0.5,
                  segment.alpha = 0.4,
                  segment.colour = "#cccccc",
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_few()

###trying out a different package

mew1 <- (assay(vsd)[select,])
colnames(mew1) <- paste(vsd$group, vsd$specification, sep = " - " )
res.pca <- prcomp(t(mew), scale = TRUE)
fviz_eig(res.pca)
uh <- rep(c(16),22)
fviz_pca_ind(PCA, geom = "point", col.ind = vsd$group, 
             pointsize = 2.5, invisible="quali") + 
  scale_shape_manual(name = "Group", values = uh) +
  scale_color_brewer(name = "Group", palette = "Set1") +
  ggtitle("") +
  theme(text = element_text(size = 4))
fviz_pca_var(PCA,
             select.var = list(names= top_genes),
             labelsize=3,
             alpha.var = 0.1,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
?fviz_pca_biplot
pal <- wes_palette("Zissou1", 100, type = "continuous")
three <- fviz_pca_biplot( PCA,
                 geom.ind = "point",
                 fill.ind = vsd$group,
                 col.ind =  "lightgrey",
                 pointsize = 2, pointshape=21,
                 palette = colscale,
                 addEllipses = FALSE,
                 # col.ind = vsd$group,
                 #variables
                 alpha.var = 0.1,
                 col.var = "contrib",
                 select.var = list(names= top_genes),
                 gradient.cols = pal,
                 #c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,
                 labelsize = 3, 
                 axes.linetype=NA)  + ggtitle("") +
  labs(x= 'PC1', y= 'PC2')+
  guides(fill=FALSE)+
  theme_few()
four <- fviz_pca_biplot( PCA,
                         axes = c(1, 4),
                          geom.ind = "point",
                          fill.ind = vsd$group,
                          col.ind =  "lightgrey",
                          pointsize = 2, pointshape=21,
                          palette = colscale,
                          addEllipses = FALSE,
                          # col.ind = vsd$group,
                          #variables
                          alpha.var = 0.1,
                          col.var = "contrib",
                          select.var = list(names= top_genes4),
                          gradient.cols = pal,
                          #c("#00AFBB", "#E7B800", "#FC4E07"),
                          repel = TRUE,
                          labelsize = 3, 
                          axes.linetype=NA)  + ggtitle("") +
  labs(x= 'PC1', y= 'PC4')+
  guides(fill=FALSE)+
  theme_few()

classic_theme <- theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())

ggplot(X) +
  geom_point(aes(-x, 0, color = vsd$group)) + 
  classic_theme + theme(axis.line.y =element_blank()) +
  scale_color_manual(values = colscale)+
  guides(color=FALSE)+
  theme_void()
?fviz_pca_biplot

##### Saving the PCA plots

p <- (one + two) / three /four 
tiff("pca2.tiff",units="in", width=13,height=17, res=300, compression = 'lzw') 
p[[1]] <- p[[1]] + plot_layout(tag_level = 'new',
                                   guides = "collect") 
p[[2]] <- p[[2]] + plot_layout(tag_level = 'new',
                               guides = "collect") 
p[[3]] <- p[[3]] + plot_layout(tag_level = 'new',
                               guides = "collect")
p +
  plot_annotation(tag_levels = c('A','1'),
                  tag_sep = '.')


dev.off()

tiff("pca2.tiff",units="in", width=9,height=4, res=300, compression = 'lzw') 
newPatch <- screePlot + ordOne + ordTwo

newPatch + 
  plot_annotation(tag_levels = 'A') &
  guides(fill=FALSE)

dev.off()
############# HeatMap

top_genesHM <- PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "1" |  PC == "2") %>%
  # for each PC
  group_by(PC) %>%
  # sort descending value
  dplyr::arrange(desc(abs(value))) %>%
  # take top 5 rows of each PC group
  dplyr::slice(1:2000) %>%
  # extract the column (gene name) from the table
  pull(column) %>%
  # retain unique gene names only
  unique()

colors <- colorRampPalette( rev(brewer.pal(9, "RdGy")) )(255)


library("factoextra")
library("magrittr")
library("cluster")
# Pairwise correlation between samples (columns)
cols.cor <- cor((assay(vsd)[top_genesHM,]), 
                use = "pairwise.complete.obs", 
                method = "pearson")

# Pairwise correlation between rows (genes)
rows.cor <- cor(t(assay(vsd)[top_genesHM,]), 
                  use = "pairwise.complete.obs", 
                  method = "pearson")

mew <- (assay(vsd)[top_genesHM,])
colnames(mew) <- paste(vsd$group, vsd$specification, sep = " - " )
#rownames(mew) <- NULL
# Plot the heatmap


# Pairwise correlation between samples (columns)
forhc <- (assay(vsd)[top_genesHM,])
colnames(forhc) <- paste(vsd$group, vsd$specification, sep = " - " )
cols.cor1 <- cor((forhc), 
                use = "pairwise.complete.obs", 
                method = "pearson")

library(dendsort)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
hc.complete = hclust(dist(t(mew)), method="complete")
hc.complete$label <- paste(vsd$group, vsd$specification, sep = " - " )
hc.complete.rows = hclust(dist((mew)), method="complete")


find_coordinates = function(n, gaps, m = 1:n){
  if(length(gaps) == 0){
    return(list(coord = unit(m / n, "npc"), size = unit(1 / n, "npc") ))
  }
  
  if(max(gaps) > n){
    stop("Gaps do not match with matrix size")
  }
  
  size = (1 / n) * (unit(1, "npc") - length(gaps) * unit("4", "bigpts"))
  
  gaps2 = apply(sapply(gaps, function(gap, x){x > gap}, m), 1, sum) 
  coord = m * size + (gaps2 * unit("4", "bigpts"))
  
  return(list(coord = coord, size = size))
}

draw_colnames = function(coln, gaps, vjust_col, hjust_col, angle_col, ...){
  coord = find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = vjust_col, hjust = hjust_col, rot = angle_col, gp = gpar(...))
  
  return(res)
}

assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames",
  ns = asNamespace("pheatmap")
)

tiff("heatMapThesis.tiff",units="in", width=10,height=10, res=300, compression = 'lzw') 
out <- pheatmap(
  mew,
  #scale = "row",
  cluster_cols=hc.complete,
  cluster_rows =hc.complete.rows,
  clustering_distance_cols = as.dist(1 - cols.cor),
  clustering_distance_rows = as.dist(1 - rows.cor),
  cutree_cols = 4,
  cutree_rows = 6,
  annotation_colors = heatmapColScale,
 #annotation_col = colAn,
 #annotation_row = rowAn,
  show_rownames = FALSE,
 #show_colnames = FALSE,
 labels_col = hc.complete$label,
  treeheight_row = 0,
 angle_col = 315,
 fontsize_col = 7,
 col=inferno(100)
#  col = inferno(length(breaksList)),
 # breaks = breaksList
)

dev.off()

################
filtered <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  #dplyr::filter(gene %in% top_genes) %>%
  select(gene, PC1, PC4) %>%
  dplyr::filter(PC4 > -0.01 & PC4 < 0.01 ) 

filtHeat <- forDen %>%
  as.tibble(rownames="gene")%>%
  select(gene,contains("EHM"))

filtHeat <- filtHeat %>%
   filter(gene %in% filtered$gene)

filtHeat <- as.data.frame(filtHeat)
rownames(filtHeat) <- filtHeat$gene
filtHeat <- filtHeat[,-1]

filtHeat <- filtHeat[apply(filtHeat, MARGIN = 1, FUN = function(x) sd(x) != 0),]

out <- pheatmap(na.omit(filtHeat),
        # scale="row",
         method = "complete",
         color = inferno(100),
        show_rownames = FALSE
        # cutree_rows = 
        # cutree_cols = 2
         )

?pheatmap
#making annotation 
#rowAn
out$tree_row$labels
reorder_idx <- match(first,second) # Saving indices for how to reorder `second` to match `first`

second[reorder_idx]  # Reordering the second vector to match the order of the first vector
second_reordered <- second[reorder_idx]  # Reordering and saving the output to a variable
rownames(clusDesignation[out$tree_row$labels,])
idx <- match(out$tree_row$labels,clusDesignation$one)
clusDesignation$one[idx]
df[match(target, df$name),]



out$tree_row
sort(cutree(out$tree_row, 6))
rowAn <- as.data.frame(sort(cutree(out$tree_row, 6)))

clusDesignation <- as.data.frame(sort(cutree(out$tree_row, 6)))
clusDesignation$one <- rownames(clusDesignation)
write.csv(clusDesignation, "heatmapClustersThesis.csv")

#colAn
colAn <- as.data.frame(vsd$group)
rownames(colAn) <- (colnames(mew))
colnames(colAn) <- c("Group")
colnames(rowAn) <- c("Cluster")

heatmapColScale <- list(
  Group = c(Adult_Heart = "#e15759", Fetal_Heart = "#9d983d",EHM="#90728f",CM="#f3a546"),
  Cluster = c('1'= "red",
              '2' = "yellow",
              '3' = "blue",
              '4' = "black",
              '5' = "pink",
              '6' = "green"))



hc <- hclust(as.dist(1- cols.cor1))
dd <- as.dendrogram(hc)
dd.reorder <- reorder(dd, 10:1)
plot(dd.reorder)
plot(hc)
forDen <- (assay(vsd)[top_genesHM,])
colnames(forDen) <- paste(vsd$group, vsd$specification, sep = " - " )
keepNeedingThis <- paste(vsd$group, vsd$specification, sep = " - " )
hc.plot = hclust(dist(t(forDen)), method="complete")
dendr    <- dendro_data(hc.plot, type="rectangle") # convert for ggplot
#clust    <- cutree(hc,k=2)                    # find 2 clusters
clust.df <- data.frame(label=keepNeedingThis, cluster=factor(vsd$group))
# dendr[["labels"]] has the labels, merge with clust.df based on label column
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
# plot the dendrogram; note use of color=cluster in geom_text(...)
tiff("Dendo-thesis.tiff",units="in", width=8,height=7, res=300, compression = 'lzw')
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
            size=3) +
  coord_flip()+
  labs(x = " ", y=" ")+
  scale_y_reverse(expand=c(0.2, 0)) + 
  scale_color_manual(values=colscale)+
  theme(axis.line.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank()
        ) +
  guides(color=FALSE)


dev.off()
tiff("Corr-thesis.tiff",units="in", width=7,height=7, res=300, compression = 'lzw') 
corrplot(cols.cor1, 
         method = "color",
         tl.col='grey30',
         addrect = 3,
         order = 'hclust',
         tl.cex = 0.5,
         number.cex = .7,
         col=inferno(200),
         #type="upper",
       #  bg = "black",
         diag = FALSE) 

dev.off()


#########################
forCorr <- (assay(vsd)[top_genesHM,])
colnames(forCorr) <- vsd$group
r <- cor(forCorr, method="pearson")
round(r,2)

s <- r[rownames(r)=="Adult_Heart",]
t <- as.data.frame(colMeans(s))
t$group <- colnames(s)
t$group <- as.factor(t$group)
colnames(t) <-c("corr", "group")

u <- t[t$group!= "Adult_Heart",]

# Using median
tiff("PearsonCorr-thesis.tiff",units="in", width=4,height=4, res=300, compression = 'lzw') 
u %>%
  mutate(class = fct_reorder(group, corr, .fun='median')) %>%
  ggplot( aes(x=reorder(group, corr), y=corr)) + 
  geom_boxplot(aes(fill = stage(group, after_scale = alpha(fill, 0.9))))+
 # geom_boxplot(fill=after_stat(),alpha=0.4) +
  ylab("Pearson Corelation\n") +
  theme(legend.position="none") +
  scale_fill_manual(values = colscale)+
  #scale_fill_tableau(palette =  "Superfishel Stone") +
  xlab("")+
  theme_few()+
  guides(fill=FALSE)

dev.off()

#######
dds = estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))

# Create tibbles including row names
metaDat <- metdata %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

dat <- normalized_counts %>% 
  as_tibble(rownames = "gene") %>%
  dplyr::filter(rownames(normalized_counts) %in% fiftyGenes)



gathered <- dat %>%
  #as.tibble(rownames = "gene") %>%
  tidyr::gather(colnames(dat)[2:77], key = "samplename", 
         value = "normalized_counts")
joined <- inner_join(metaDat,gathered )

joined <- joined %>%
  filter(group != "ipsc",
         group != "Fib",
         group != "Rh",
         fromWhere != "Cyg")


xb <- joined %>%
  group_by(group, gene) %>%
  dplyr::summarize(median = median(normalized_counts, na.rm=TRUE),
                   sd = sd(normalized_counts, na.rm=TRUE),
                   mean = mean(normalized_counts, na.rm=TRUE))

xb <- inner_join(xb,joined)

lb <- xb %>%
  filter(group == "Adult_Heart") %>%
  mutate(reOrderedGenes = reorder(gene,median)) %>% 
  select(group,gene, reOrderedGenes, sd)

xb <- xb %>% inner_join(lb, by="gene" )


xb$group <- as.factor(xb$group.x)

xb$group <- fct_relevel(xb$group,"Adult_Heart",
                         "Fetal_Heart",
                         "EHM",
                         "CM")


upfiftyGenes <- ggplot() +
  geom_point(data=xb, aes(x = reOrderedGenes, y = normalized_counts, color=group),alpha=0.1,  size=1.5) +
  geom_errorbar(data=xb, mapping = aes(x = reOrderedGenes, y = median,  ymin = mean - sd.x, ymax = mean + sd.x, color=group),size=0.2, width=.2, alpha=0.3)+
  geom_point(data=xb, aes(x=reOrderedGenes, y=median, color= group), size=4)+
#  geom_line(data=xb, aes(x=reOrderedGenes, y=median, color= group, group=group))+
  scale_y_log10()+
  labs(x=" ",
       y=" ",
       title = "General CM Marker Genes") +
  #scale_y_log10(label = ff_denom())+
  
  #ylab("log10 Normalized Counts") +
  #ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  scale_colour_manual(values=colscale)+
  scale_fill_manual(values=colscale)+
  guides(color=FALSE,
         fill=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        legend.text=element_text(size=12),
        panel.grid= element_blank(),
        legend.title = element_blank())

tiff("GeneExp-thesis01.tiff",units="in", width=17,height=12, res=300, compression = 'lzw') 
onepat <- (fg + sg) +
  plot_layout(widths = c(1,2))
onepat / oxre / upfiftyGenes +
  plot_layout(guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  ylab("Normalized Counts") &
  theme(plot.title=element_text(size=16),
        legend.text = element_text(size=18),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=13),
        axis.title.y = element_text(size=12))
dev.off()
# Plot
library(GGally)
dat <- normalized_counts %>% 
    as_tibble(rownames = "gene") %>%
    dplyr::filter(rownames(normalized_counts) %in% gen)
  
tryout <- as.data.frame(sort(cutree(out$tree_row, 3)))
tryout$gene <- rownames(tryout)
colnames(dat)
colnames(tryout) <- c("cluster","gene")
trieedout <- inner_join(tryout, dat)
c
gathered <- dat %>%
  #as.tibble(rownames = "gene") %>%
  tidyr::gather(colnames(dat)[2:77], key = "samplename", 
                value = "normalized_counts")

joined <- inner_join(metaDat,gathered )
joined <- joined %>%
  filter(group != "ipsc",
         group != "Fib",
         group != "Rh",
         fromWhere != "Cyg")

joined$cluster <- as.character(joined$cluster)
joined$group <- fct_relevel(joined$group,
                            "CM",
                            "Fetal_Heart",
                            "EHM",
                            "Adult_Heart",
                        )

letsSee<- joined %>%
  group_by(group,gene) %>%
  summarize(mean = mean(normalized_counts, na.rm=TRUE)) %>%
 # mutate_if(is.numeric, scale) %>%
  spread(key="group", value="mean")

letsSee[2:5]<- log10(letsSee[2:5]+1)

ggparcoord(letsSee,
             columns = 2:5, groupColumn = 1, alphaLines = 0.3
  )+
  #scale_color_npg()+
  #scale_y_continuous(labels = ff_denom(), expand = c(0, 0),
   #                  limits = c(0, max(letsSee$Fetal_Heart) * 1.05) +
  guides(color=FALSE)+
  theme_few()

library(plotly)
ggplotly(luh)

# Data set is provided by R natively
data <- iris

# Plot
ggparcoord(data,
           columns = 1:4, groupColumn = 5
)   

top_genes

?sortLvlsByVar.fnc
onepat <- (upHH |upSC ) +
  plot_layout(widths = c(1,2))






dev.off()


perc <- c("75","81","93","83")
grp <- c("a","b","c","d")
dat <- data.frame(grp,perc)
dat$perc <- as.numeric(dat$perc)
res.aov <- aov(perc~grp, data=dat)
summary(res.aov)
