#now read in the writted table and the metadata of the same 
raw_counts_df <- read.table("data/rawCounts_onlyCM_thesis.txt", sep ="\t")


coldat <- read.table("data/onlyCM_coldat.txt", sep ="\t")
rownames(coldat) <- coldat[,1]


hist(raw_counts_df, notch=TRUE)
coldat$Specification <- as.factor(coldat$Specification)

coldat$Folder <- as.factor(coldat$Folder)

coldat$Group <- as.factor(coldat$Group)

dds <- DESeqDataSetFromMatrix(countData = raw_counts_df, colData = coldat, design=~Folder)
dds <- dds[ rowSums(counts(dds)) > 5, ]

dds = estimateSizeFactors(dds)
sizeFactors(dds)
mean(raw_counts_df$p722s3C190604_Tiburcy_S37_L005_R1_001.bam)
mean(raw_counts_df$p637sDiff6CM_Zimmermann_S46_L005_R1_001.bam)
wel <- raw_counts_df$p637sDiff3CM_Zimmermann_S44_L005_R1_001.bam == "0"
which(wel == "TRUE")
table(wel)["TRUE"]
length(which(wel))
vsd <- vst(dds)
boxplot(assay(rl), notch=TRUE)
rl <- rlog(dds)
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
#Z-scale normalization
logcounts <- t( apply(logcounts, 1 , function(x) scale( x , center = T, scale = T)))
boxplot(logcounts, notch = TRUE)
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

one <- ggplot(x, aes(x$PC1, x$PC2, color = p$metadata$Specification, shape = p$metadata$Group))+
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
one
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

PC1 <- as.data.frame(p$loadings[,1], drop=FALSE)
PC1$gene <- row.names(PC1)
PC1 <- as.data.frame(abs(PC1), drop=FALSE)
a <- get_top_Geneloadings(p$loadings,1,1:100)
b <- get_top_Geneloadings(p$loadings,2,1:10)
c <- get_top_Geneloadings(p$loadings,3,1:20)
meh <- c(a,b,c)
a
gene_loadings <- p$loadings %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::filter(gene %in% a) 
gene_loadings <- gene_loadings[,1]
gene_loadings[order("PC1"),]
two <- ggplot(gene_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_label_repel(data = gene_loadings, 
                   aes(x = PC1, y = PC2, label = gene))
dev.off()
tiff("99_01_onlyCM_LP_PCA.tiff",units="in", width=15,height=8, res=300, compression = 'lzw') 
one | two
dev.off()
screeplot(p)
findElbowPoint(p$variance)

tiff("99_01_sub_pairsplot.tiff",units="in", width=10,height=8, res=300, compression = 'lzw') 
pairsplot(p,
          components = getComponents(p, c(1:3)),
          triangle=FALSE, trianglelabSize = 12,
          hline=0, vline = 0,
          pointSize = 3,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby='Specification')


dev.off()
plotloadings(p,
             components = getComponents(p, c(1,3)),
             rangeRetain = 0.0001)

#######Differential expression b/w the 4 CMs 

dds <- DESeq(dds)
dds$Specification
res_tableOE_unshrunken <- results(dds, alpha = 0.05)
summary(res_tableOE_unshrunken)
#res_tableOE <- lfcShrink(dds,,res=res_tableOE_unshrunken)

res_tableOE_unshrunken %>% data.frame() %>% View()
### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
res_tableOE_unshrunken <- res_tableOE_unshrunken %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sigOE <- res_tableOE_unshrunken %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
plotCounts(dds, gene="TNNI3", intgroup="Specification", returnData = TRUE) 

########## HEATMAP of relevant genes
genes <- c("MYH7","MYBPC3","ACTN2","TNNI3","MYL2","ATP2A2","MYH6",
           "TNNI1","MYL7","NKX2-5","GATA4")
new <- logcounts[rownames(logcounts) %in% genes,]
colnames(logcounts) <- dds$Specification

tiff("99_heatmap_onlyCM.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
pheatmap(na.omit(new),
         scale="row",
         color=viridis(20),
         #main = "",
         #cutree_cols = 7,
         cluster_rows=TRUE,
         cluster_cols = TRUE,
         #color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         show_rownames = TRUE,
         show_colnames = TRUE)
dev.off()
fig3d <- c("FOXC1","FN1","ISL1","TAL1","CDH5",
           "KDR","FGF10","BMP4","THY1","PDGFRB",
           "HOXA3","HOXA1","PITX2","TBX18","GJA1",
           "TNNI1","MYL7","MYH6","NKX2-5","IRX4",
           "HCN4","FOXC2")
new <- logcounts[rownames(logcounts) %in% fa,]
tiff("99_02heatmap-fa_onlyCM.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
pheatmap(na.omit(new),
         scale="row",
         color=viridis(20),
         #main = "",
         #cutree_cols = 7,
         cluster_rows=TRUE,
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
###########
colnames(logcounts) <- dds$Specification
new <- logcounts[select,]
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
         p.mat = res2$P, sig.level = 0.005, insig = "label_sig")
dev.off()
# Insignificant correlations are leaved blank
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")
#install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")

chart.Correlation(logcounts, histogram=TRUE, pch=19)




