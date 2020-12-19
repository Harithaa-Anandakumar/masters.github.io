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
library(EDASeq)
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
library(biomaRt)
library(rafalib)
library(PCAtools)
library(purrr)
library(readr)
library(dplyr)
library(readxl)
library(FactoMineR)
library(factoextra)
library(MGFR)
library(ggbeeswarm)
library(d3heatmap)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(sva)
library(gridExtra)
library(grid)
library(gghighlight)


new <- vsd[,vsd$GMP != "1"]
dds$Group <- as.factor(dds$Group)
vsd_filtered <- vsd[,vsd$GMP == "GMP"]
rlog <- rlog(dds_filtered)
assay(rlog) <- limma::removeBatchEffect(assay(rlog), dds$Folder)

sub_dds <- vsd[,vsd$Group == "CM"]
prcomp(assay(vsd_filtered))
plot_pca <- plotPCA(vsd, intgroup= c("Group","Specification"))
ggplotly(p)
metadata <- colData(vsd_filtered)
x <- (assay(vsd_filtered))
dim(x)
class(x)
x
assay(vsd)
x <- as.data.frame(x)
dev.off()
plot_pca | loadings_plot
library(tidyverse)
library(broom)


PCA$x

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
  dplyr::slice(1:20) %>%
  # extract the column (gene name) from the table
  pull(column) %>%
  # retain unique gene names only
  unique()

#######
new <- PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "1" |  PC == "2" | PC == "3") %>%
  # for each PC
  group_by(PC) %>%
  # sort descending value
  dplyr::arrange(desc(abs(value))) %>%
  # take top 5 rows of each PC group
  dplyr::slice(1:50)
class(new)
new <- PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "4") %>%
  transform(value, squared = value^2)

sum(new$squared)

sampleDists <- dist((PCA$x[,1:5]))
Dists <- dist(t(logcounts))
colnames(logcounts) <- samples
Dists
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$Group, vsd$Specification, sep = " - " )
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
tiff("dist-mat2",units="in", width=7,height=7, res=300, compression = 'lzw')
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

dev.off()


ggplot(new, aes(x = column, y = value)) +
  geom_dotplot()



top_genes
gene_loadings <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::filter(gene %in% top_genes) 

loadings_plot <- ggplot(gene_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_label_repel(data = gene_loadings, 
                   aes(x = PC1, y = PC2, label = gene))

zwei
loadings_plot
# Adjust some aspects of each plot
plot_pca <- plot_pca + 
  coord_fixed(ratio = 1.0) + 
  labs(tag = "A", title = "PC scores") + 
  theme(legend.position = "none") 


loadings_plot <- loadings_plot + 
  coord_fixed(ratio = 0.5) + 
  labs(tag = "B", title = "PC loadings")

dev.off( )
drei | loadings_plot

plot_grid(pca_plot, loadings_plot, ncol = 2, align = "h")
#using PCAtools
p <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)
screeplot(p)
pca_plot <- biplot(p)
plotloadings(p)
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 3.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)

pairsplot(p)
dev.off()
#looking at loading scores to determine which genes have the largest effect on where samples are plotted in PCA plot
loading_scores <- PCA$rotation[,2] #prcomp func calles loading scores rotation, there are LS for each principal component
#abs() gives both sets of genes
rownames(loading_scores) <- rownames(p$loadings)
# genes pushing samples to left - large negative value and to the rt - large pos value
gene_scores <- abs(loading_scores)

#now sorting magnitude from high to low
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
(gene_score_ranked)
#top 10 genes
top_20_genes <- c(top_10_genes,meh)

meh<- names(gene_score_ranked[1:10])

#looking at loading scores to determine which genes have the largest effect on where samples are plotted in PCA plot
loading_scores <- PCA$rotation[,1] #prcomp func calles loading scores rotation, there are LS for each principal component
#abs() gives both sets of genes
rownames(loading_scores) <- rownames(p$loadings)
# genes pushing samples to left - large negative value and to the rt - large pos value
gene_scores <- abs(loading_scores)

#now sorting magnitude from high to low
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
(gene_score_ranked)
#top 10 genes
top_10_genes <- names(gene_score_ranked[1:10])
top_10_genes
write.csv(top_10_genes,"")
write.table(top_10_genes, file = "mtcars.txt", sep = " ",
            row.names = FALSE, quote = FALSE)



#########################################


sub_vsd <- vsd[,vsd$GMP !="non"]
plotPCA(sub_vsd, intgroup = "GMP")
plotPCA(vsd, intgroup = "Folder")
dds$GMP <- factor(dds$GMP, levels = c("GMP","RH"))
dds$GMP <- droplevels(dds$GMP)
dds$GMP <- as.factor(dds$GMP)
new <- DESeq(dds)
res <- results(new)
res
?results
resultsNames(new)
# generate results table for GMP vs RH
res <- results(new, name="GMP_RH_vs_GMP")
resLFC <- lfcShrink(new, coef="GMP_RH_vs_GMP", type="apeglm")
# because we are interested in treated vs untreated, we set 'coef=2'
resNorm <- lfcShrink(new, coef=2, type="normal")
resAsh <- lfcShrink(new, coef=2, type="ashr")
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm = TRUE)
resSig <- subset(resOrdered, padj < 0.005)
plotMA(resLFC, ylim=c(-2,2))
write.csv(resSig,"GMP_vs_RH.csv")

## 

p <- pca(na.omit(assay(vsd_filtered)), metadata=(colData(vsd_filtered)),removeVar = 0.1)
p <- pca(na.omit(assay(vsd)), metadata=(colData(vsd)),removeVar = 0.1)
p <- pca(na.omit(assay(dds)), metadata=(colData(dds)),removeVar = 0.1)
ggplot(p$rotated, aes(x =p$rotated[,"PC1"], y = dds$Folder, color = dds$Specification)) +
  geom_quasirandom(groupOnX = TRUE) +
  geom_point(size=4)+
  theme_bw()+
  theme(legend.position = "none") +
  xlab("PC1") + ylab("Cell Types -clubbed") +
  ggtitle("Cells ordered by principal components")+
  scale_fill_pander()
