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
f_files<- list.files("/Users/bhuvaneswari/Desktop/thesis/self_contained_thesis/data",
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

#write.table(raw_counts_df,"rawCounts_thesis.txt", sep="\t", quote = FALSE)

#now read in the writted table and the metadata of the same 
raw_counts_df <- read.table("data/rawCounts_thesis.txt", sep ="\t")

#remove first column of geneID 

raw_counts_df <- as.data.frame(raw_counts_df)
raw_counts_df <- raw_counts_df[,-1]

coldat <- read.csv("data/coldat_thesis.csv")
three <- read_xlsx("../GeneList.xlsx")
three <- three$`Table S6:  Genes expressed at greater than 1 FPKM in 75% of cardiomyocytes at P0  used to produce Figure 4B,  5A, 5B, 6A`
three <- toupper(three)
#create a DESEQ object 
dds <- DESeqDataSetFromMatrix(countData = raw_counts_df, colData = coldat, design=~ meh)
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
#DON"T CORRECT _ IT"S WRONG

#convert to log counts with an addition of one
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
colnames(logcounts) <- dds$Specification
#Z-scale normalization
logcounts <- t( apply(logcounts, 1 , function(x) scale( x , center = T, scale = T)))
boxplot(logcounts, notch = TRUE)
dev.off()
# perform pca using prcomp func, choosing the top n number of genes 
ntop=2000
Pvars <- rowVars(logcounts)
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]

##creating PCA object based on Kevin B's PCAtools.
colnames(logcounts) <- rownames(colData(dds))
q <- pca(na.omit(logcounts), metadata=(colData(dds)),removeVar = 0.1)
x <- p$rotated
rownames(x) <- NULL
comp <- p$loadings

#percent variation each component accounts for
percentVarp <- round(100*p$sdev^2/sum(p$sdev^2),1)

one <- ggplot(x, aes(x$PC1, x$PC2, color = p$metadata$Group, shape = p$metadata$Group))+
  geom_point(size =3) +
  #ggtitle("PCA-using all genes") +
  theme(plot.title = element_text(hjust=0, size=8))+
  geom_label_repel(data = x, 
                  aes(PC1, PC2),
                  label = p$metadata$Specification) +
  xlab(paste0("PC1: ", round(percentVarp[1],digits=2),"% variance"))+
  ylab(paste0("PC2: ", round(percentVarp[2],digits = 2), "% variance"))+
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

a <- get_top_Geneloadings(p$loadings,1,1:20)
b <- get_top_Geneloadings(p$loadings,2,1:20)
c <- get_top_Geneloadings(p$loadings,3,1:20)
meh <- c(a,b,c)

gene_loadings <- p$loadings %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::filter(gene %in% meh) 

two <- ggplot(gene_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_label_repel(data = gene_loadings, 
                   aes(x = PC1, y = PC2, label = gene))
tiff("99_01full_LP_PCA.tiff",units="in", width=15,height=8, res=300, compression = 'lzw') 
one | two
dev.off()
#######


sarcomeric_assembly <- c("MYH6", "CAPN3", "NEB", "EDN","NKX2-5")

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
dev.off()
tiff("99_full_LP_PCA.tiff",units="in", width=15,height=8, res=300, compression = 'lzw') 
one | loadings_plot 
dev.off()



########heat map !
genes <- c("MYH7","MYBPC3","ACTN2","TNNI3","MYL2","ATP2A2","MYH6",
           "TNNI1","MYL7","NKX2-5","GATA4")
new <- logcounts[rownames(logcounts) %in% three,]
colnames(new) <- dds$Specification
onlyCM <- dds[,dds$Group =="CM"]
onlyCM <- vst(onlyCM)
new <- assay(onlyCM)[rownames(assay(onlyCM)) %in% three,]
colnames(new) <- onlyCM$Specification
tiff("99_heatmap_onlyCM.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
out <- pheatmap(na.omit(new),
         scale="row",
         #color=viridis(20),
         #main = "",
         #cutree_cols = 7,
         cluster_rows=TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
         cex=1, clustering_distance_rows="correlation", 
         cex=1,
         clustering_distance_cols="correlation", 
         clustering_method="complete", border_color=FALSE,
         show_rownames = FALSE,
         show_colnames = TRUE)

out
out <- pheatmap(new, 
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
         cex=1, clustering_distance_rows="correlation", cex=1,
         clustering_distance_cols="correlation", clustering_method="complete", border_color=FALSE)

summary(out)

sort(cutree(out$tree_row, k=5))
plot(out$tree_row)
abline(h=7, col="red", lty=2, lwd=2)
clusDesignation <- cutree(as.hclust(out$tree_row), 5)
names(clusDesignation[clusDesignation==1])
names(clusDesignation[clusDesignation==2])
names(clusDesignation[clusDesignation==3])
names(clusDesignation[clusDesignation==4])
names(clusDesignation[clusDesignation==5])
names(clusDesignation[clusDesignation==6])
names(clusDesignation[clusDesignation==7])

dev.off()

fig3d <- c("FOXC1","FN1","ISL1","TAL1","CDH5",
           "KDR","FGF10","BMP4","THY1","PDGFRB",
           "HOXA3","HOXA1","PITX2","TBX18","GJA1",
           "TNNI1","MYL7","MYH6","NKX2-5","IRX4",
           "HCN4","FOXC2")
new <- logcounts[rownames(logcounts) %in% glycolysis,]
tiff("99_01heatmap-glycolysis_onlyCM.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
pheatmap(na.omit(new),
         scale="row",
         color=viridis(20),
         #main = "",
         #cutree_cols = 7,
         cluster_rows=TRUE,
         clustering_distance_rows="euclidean",
         cex=1,
         clustering_distance_cols="euclidean", 
         clustering_method="complete", 
         border_color=FALSE,
         cluster_cols = TRUE,
         #color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         show_rownames = TRUE,
         show_colnames = TRUE)
dev.off()
fa <- c("ACADM","CPT2","CPT1B","ETFDH","HADHA","HADHB","MAPK14")
oxre <- c("HIBADH","IMPDH2","NDUFA10","NDUFA11","NDUFA12","NDUFB6","NDUFS1","NDUFV1","ACAD9","ACADM","ADHFE1","ALDH4A1","ALDH5A1",
          "ALDH6A1","ALDH7A1","AASS","AJFM1","CAT","CRYZ","CYP27A1","CYBSA","CYB5A","COX681","DHTKD1","DHFR","DLD","ETFA","ETFB","ETFDH",
          "FTH1","GLUD1","GCDH","GPD2","GRHPR","HADH","HSD17B10","IDH2","MDH2","OGDH","PRDX3",
          "PCYOX1","PYROXD2","PDHA1","PDHB","SLC25A12","SLC25A13","SUOX","TMLHE","UQCRF51","UQCRFS1")
glycolysis <- c("BPGM","TPI1","ALDOA","ALDOC","ENO1","ENO2","GPI","GAPDH","HK2","LDHA","PFKL","PFKP","PGK1","PGAM1")





#####don't run stuff below -- issue 
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

###########
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
distancem <- dist(top_genes_exprs.mx)
hclust_completem <- hclust(distancem, method = "complete")
dendcompletem <- as.dendrogram(hclust_completem)
heatmap(new, Rowv=dendcompletem, Colv=NA, scale="column")
############################
# Calculate correlations between samples.
colnames(top_genes_exprs.mx) <- dds$Specification
cr = cor(top_genes_exprs.mx, method='pearson')

# Calculate correlations between genes.
gene_cor = cor(t(top_genes_exprs.mx), method='pearson')

# Find gene distances:
gene_dist = as.dist(1-gene_cor)

# hclust the data
hc_genes = hclust(gene_dist, method='complete')
hc_samples = hclust(as.dist(1-cr), method='complete')

# heat map colors and breaks.
myheatcol = colorRampPalette(c('red', 'yellow'))(n = 75)
quantBrks = quantile(top_genes_exprs.mx, c(0.03, 0.97))

# Make heatmap
heatmap.2(top_genes_exprs.mx,
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
##########################################################
dds$Group

dds <- DESeq(dds)
res <- results(dds, contrast = c('meh','Ventricular CM','CM-GMP'))
head(res1)
summary(res1)
# Remove genes with no p-value
uninformative_genes <- is.na(res$padj) 
result.df <- res[!uninformative_genes,]
dds$meh
res <- lfcShrink(dds,
                  contrast = c('meh','Ventricular CM','CM-GMP'), res=res, type ="apeglm")
EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 8))

# Set cut-offs for top genes
log_fc_cutoff <- 2
adj_p_cutoff <- 0.05
result.df <- as.data.frame(result.df)
result.df$gene <- row.names(result.df)
# Extract results for top genes
top_genes.df <- result.df %>%
  dplyr::filter(abs(log2FoldChange) > log_fc_cutoff, padj < adj_p_cutoff)


write.csv(top_genes.df, "DE_RH.csv",sep="\t", quote = FALSE)
# Get expression values for the top genes
top_genes_exprs.mx <- counts(dds, normalized=TRUE)[top_genes.df$gene, ] # Change gene IDs to gene names

# Plot heatmap
num_genes <- nrow(top_genes.df)
colnames(top_genes_exprs.mx) <- dds$Specification
main <- paste(num_genes, "genes with log2(FC) >", log_fc_cutoff, "and adj.P <", adj_p_cutoff,  "\n")
par(oma=c(0,0,3,0)) # make space for the main header
heatmap(top_genes_exprs.mx, main=main)

main  

library(clusterProfiler)
geneList <- top_genes.df$log2FoldChange
names(geneList) <- as.character(top_genes.df$gene)
geneList <- sort(geneList, decreasing = TRUE)
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
merge(geneList, gene.df)
gene <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)
library(org.Hs.eg.db)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

edo2 <- gseNCG(geneList, nPerm=10000)

#######################

