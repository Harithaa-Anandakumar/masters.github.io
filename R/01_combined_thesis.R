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
library(gridExtra)
library(grid)
library(gghighlight)
library(plotly)


f_files<- list.files("counts",
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
raw_counts_df <- raw_counts_df[,-1]

sample <- colnames(raw_counts_df)
write.csv(sample,"combined_colDat_thesis.csv")

coldat <- read.csv("combined_colDat_thesis.csv")
rownames(coldat) <- coldat$x
three <- read_xlsx("../GeneList.xlsx")
three <- three$`Table S6:  Genes expressed at greater than 1 FPKM in 75% of cardiomyocytes at P0  used to produce Figure 4B,  5A, 5B, 6A`
three <- toupper(three)

write.table(raw_counts_df,"combined_rawCounts_thesis.txt", sept="\t",quote=F)

dds <- DESeqDataSetFromMatrix(countData = raw_counts_df, colData = coldat, design=~ Folder)
dds <- dds[ rowSums(counts(dds)) > 5, ]
vsd <- vst(dds)
p <- plotPCA(vsd, intgroup="meh")
ggplotly(p)
#The following code estimates size factors to account for differences in sequencing depth. 
#So the value are typically centered around 1. 
#If all the samples have exactly the same sequencing depth, you expect these numbers to be near 1. 
dds = estimateSizeFactors(dds)
sizeFactors(dds)
colSums(counts(dds))
#Plot column sums according to size factor
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

#convert to log counts with an addition of one
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
colnames(logcounts) <- dds$Specification
#Z-scale normalization
logcounts <- t( apply(logcounts, 1 , function(x) scale( x , center = T, scale = T)))
boxplot(logcounts, notch = TRUE)
boxplot(assay(vsd), notch=TRUE)
dev.off()
# perform pca using prcomp func, choosing the top n number of genes 
ntop=1000
Pvars <- rowVars(logcounts)
selectn <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]

##creating PCA object based on Kevin B's PCAtools.
colnames(logcounts) <- rownames(colData(dds))
p <- pca(na.omit(logcounts), metadata=(colData(dds)),removeVar = 0.1)
x <- p$rotated
rownames(x) <- NULL
comp <- p$loadings

#percent variation each component accounts for
percentVarp <- round(100*p$sdev^2/sum(p$sdev^2),1)

one <- ggplot(x, aes(x$PC1, x$PC2, color = p$metadata$meh, shape = p$metadata$Group))+
  geom_point(size =3) +
  #ggtitle("PCA-using all genes") +
  theme(plot.title = element_text(hjust=0, size=8))+
  geom_label_repel(data = x, 
                   aes(PC1, PC2),
                   label = p$metadata$meh) +
  xlab(paste0("PC1: ", round(percentVarp[1],digits=2),"% variance"))+
  ylab(paste0("PC2: ", round(percentVarp[2],digits = 2), "% variance"))+
  guides(color=FALSE, shape=FALSE)+
  theme_few()
one
new <- logcounts[selectn,]
colnames(new) <- dds$meh
onlyCM <- dds[,dds$Group =="CM"]
onlyCM <- vst(onlyCM)
new <- assay(onlyCM)[rownames(assay(onlyCM)) %in% three,]
colnames(new) <- onlyCM$Folder
tiff("99_heatmap_onlyCM.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
assay(onlyCM) <- limma::removeBatchEffect(assay(onlyCM), onlyCM$Folder) 
pheatmap(na.omit(new),
                scale="row",
                #color=viridis(20),
                #main = "",
                #cutree_cols = 7,
                cluster_rows=TRUE,
                cluster_cols = TRUE,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                 clustering_distance_rows="correlation", 
                clustering_distance_cols="correlation", 
                clustering_method="ward", border_color=FALSE,
                show_rownames = FALSE,
                show_colnames = TRUE)

install.packages("corrr")
library(corrr)
t(new) %>% correlate() %>% network_plot(min_cor=0.6)


colnames(logcounts) <- dds$Specification
new <- logcounts[select,]
new <- top_genes_exprs.mx
res <- cor(new)
library(corrplot)
corrplot(res, order="hclust", col = brewer.pal(n = 8, name = "RdBu"),
         tl.col='black')
# Get some colors
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE)

###
library(Hmisc)
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
res2 <- rcorr(as.matrix(new))
flattenCorrMatrix(res2$r, res2$P)
# Insignificant correlation are crossed
corrplot(res2$r, type="upper", order="FPC", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")
dev.off()
# Insignificant correlations are leaved blank
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")
#install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")

chart.Correlation(new, histogram=TRUE, pch=19)



#########
distancem <- dist(new)
hclust_completem <- hclust(distancem, method = "complete")
dendcompletem <- as.dendrogram(hclust_completem)
heatmap(new, Rowv=dendcompletem, Colv=NA, scale="column")
############################
# Calculate correlations between samples.
colnames(new) <- dds$new
cr = cor(new, method='pearson')

# Calculate correlations between genes.
gene_cor = cor(t(new), method='pearson')

# Find gene distances:
gene_dist = as.dist(1-gene_cor)

# hclust the data
hc_genes = hclust(gene_dist, method='complete')
hc_samples = hclust(as.dist(1-cr), method='complete')

# heat map colors and breaks.
myheatcol = colorRampPalette(c('red', 'yellow'))(n = 75)
quantBrks = quantile(new, c(0.03, 0.97))

# Make heatmap
heatmap.2(new,
          key.xlab="log2 centered TMM-norm CPM",
          margins=c(10,5),
          lmat=rbind(c(4,3),c(2,1)),
          lhei=c(1.3,5),
          lwid=c(2.5,5),
          labRow=NA,
          dendrogram='both',
          Rowv=as.dendrogram(hc_genes),
          Colv=as.dendrogram(hc_samples),
          col=myheatcol,
          scale="row",
          density.info="none",
          trace="none",
          key=TRUE,
          keysize=1.2,
          cexCol=1,
          breaks=seq(quantBrks[1], quantBrks[2], length=76))
dev.off()
library(GGally)
ggcorr(new, nbreaks=8, palette='RdGy', label=TRUE, label_size=5, label_color='white')
p <- plotPCA(onlyCM, intgroup="new")
ggplotly(p)
