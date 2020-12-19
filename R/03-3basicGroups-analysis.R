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


#import all .txt files from  a certain folder -- these are rawCount files produced from 
#some abundance estimation software such as FeatureCounts
f_files<- list.files("/Users/bhuvaneswari/Desktop/thesis/self_contained_thesis/data/subgroup/",
                     pattern = "*.txt", full.names = T)

read_in_feature_counts<- function(file){
  cnt<- read_tsv(file, col_names = T, comment = "#")
  cnt<- cnt %>% dplyr::select(-Chr,-Start, -End, -Strand, -Length)
  return(cnt)
}
raw_counts<- map(f_files, read_in_feature_counts)
raw_counts_df<- purrr::reduce(raw_counts, inner_join) 
raw_counts_df <- as.data.frame(raw_counts_df)
row.names(raw_counts_df) <- raw_counts_df$Geneid

#perform the above only once and save the combined rawCounts file.

write.table(raw_counts_df,"rawCounts_subgroup_thesis.txt", sep="\t", quote = FALSE)
write.csv(colnames(raw_counts_df),"coldat_subgroup.csv")


#now read in the writted table and the metadata of the same 
raw_counts_df <- read.table("rawCounts_subgroup_thesis.txt", sep ="\t")

#remove first column of geneID 
raw_counts_df <- raw_counts_df[,-1]

coldat <- read.csv("/Users/bhuvaneswari/Desktop/thesis/self_contained_thesis/coldat_subgroup.csv")
rownames(coldat) <- coldat[,1]
coldat$Specification <- as.factor(coldat$Specification)
coldat$Folder <- as.factor(coldat$Folder)
coldat$Group <- as.factor(coldat$Group)
#create a DESEQ object 
dds <- DESeqDataSetFromMatrix(countData = raw_counts_df, colData = coldat, design=~Folder)
dds <- dds[ rowSums(counts(dds)) > 5, ]

#The following code estimates size factors to account for differences in sequencing depth. 
#So the value are typically centered around 1. 
#If all the samples have exactly the same sequencing depth, you expect these numbers to be near 1. 
dds = estimateSizeFactors(dds)
sizeFactors(dds)
colSums(counts(dds))
#Plot column sums according to size factor
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

#The argument normalized equals true, divides each column by its size factor.
#assay(dds) <- limma::removeBatchEffect(assay(dds), dds$Folder)

#convert to log counts with an addition of one
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )

#Z-scale normalization
logcounts <- t( apply(logcounts, 1 , function(x) scale( x , center = T, scale = T)))

# perform pca using prcomp func, choosing the top n number of genes 
ntop=2000
Pvars <- rowVars(logcounts)
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]
colnames(logcounts) <- coldat$X

p <- pca(na.omit(logcounts), metadata=(colData(dds)),removeVar = 0.1)
x <- p$rotated
rownames(x) <- NULL
comp <- p$loadings

#percent variation each component accounts for
percentVarp <- round(100*p$sdev^2/sum(p$sdev^2),1)

one <- ggplot(x, aes(x$PC1, x$PC3, color = p$metadata$Group, shape = p$metadata$Group))+
  geom_point(size =3) +
  #ggtitle("PCA-using all genes") +
  theme(plot.title = element_text(hjust=0, size=8))+
  geom_label_repel(data = x, 
                   aes(PC1, PC3),
                   label = p$metadata$Specification) +
  xlab(paste0("PC1: ", round(percentVarp[1],digits=2),"% variance"))+
  ylab(paste0("PC3: ", round(percentVarp[3],digits = 2), "% variance"))+
  guides(color=FALSE, shape=FALSE)+
  theme_few()

## Currently works only with the pca object created by PCAtools 
## loadings dataframe -- p$loadins, 
##loadings column -- a number denoting the PC (1 or 2 etc)
##range of genes -- top 20 (1:20, etc
get_top_Geneloadings <- function(loadings_dataframe, loadings_column, range_of_genes) {
  loading_scores <- (loadings_dataframe[,loadings_column, drop=FALSE])  
  loading_scores <- abs(loading_scores)
  loading_scores$gene.names <- rownames(loading_scores)
  gene_score_ranked <- loading_scores[order(loading_scores[,1]),]
  a <- (gene_score_ranked[range_of_genes,2])
  return(a)
}

a <- get_top_Geneloadings(p$loadings,1,1:10)
b <- get_top_Geneloadings(p$loadings,2,1:10)
c <- get_top_Geneloadings(p$loadings,3,1:20)
meh <- c(a,c)

gene_loadings <- p$loadings %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::filter(gene %in% meh) 

two <- ggplot(gene_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_label_repel(data = gene_loadings, 
                   aes(x = PC1, y = PC2, label = gene))
dev.off()
tiff("99_01_sub_LP_PCA.tiff",units="in", width=15,height=8, res=300, compression = 'lzw') 
one | two
dev.off()
screeplot(p)
findElbowPoint(p$variance)

tiff("99_01_sub_pairsplot.tiff",units="in", width=10,height=8, res=300, compression = 'lzw') 
pairsplot(p,
          components = getComponents(p, c(1:5)),
          triangle=FALSE, trianglelabSize = 12,
          hline=0, vline = 0,
          pointSize = 1,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby='Group')


dev.off()
plotloadings(p, 
             rangeRetain = 0.001)







###################
meh <- na.omit(logcounts)
meh <- na.omit(assay(sub))
colnames(meh) <- sub$Specification
library(corrplot)
m <- cor(t(meh)
         library(RColorBrewer)
         
         tiff("corr-plot_skeletal_muscle",units="in", width=10,height=10, res=300, compression = 'lzw') 
         corrplot(m, order="FPC", col = brewer.pal(n = 8, name = "RdBu"),
                  tl.col='black')
         dev.off()

colnames(logcounts) <- coldat$Specification
Dists <- dist(t(assay(vsd)[select,]))

sampleDistMatrix <- as.matrix(Dists )
rownames(sampleDistMatrix) <- paste( dds$Group, dds$Specification, sep = " - " )
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
tiff("dist-mat2",units="in", width=7,height=7, res=300, compression = 'lzw')
pheatmap(sampleDistMatrix,
         clustering_distance_rows = Dists,
         clustering_distance_cols = Dists,
         col = colors)

dev.off()



