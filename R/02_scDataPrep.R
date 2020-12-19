library(cellrangerRkit)
library(dplyr)
library(Matrix)
library(gdata)
library(R.utils)
library(scran)
library(limSolve)
library(reshape2)
library(ggplot2)
library(dynamicTreeCut)
library(scater)
library(edgeR)

#read in data files 
data <- read.table(file ="/Users/bhuvaneswari/Desktop/Thesis/E-MTAB-6268/Day15Rep2_ReadMappedCount.txt")

data2 <- read.table(file ="/Users/bhuvaneswari/Desktop/Thesis/E-MTAB-6268/Day30Rep2_ReadMappedCount.txt",
                    sep="\t")
object.size(data)

data <- as.matrix(data)
matSparse <- as(data, "sparseMatrix")
object.size(matSparse)

sce_HiPSC <- SingleCellExperiment(list( counts =  matSparse))

# Find mitochondrial genes in our HiPSC dataset
is.mito_HiPSC <- grepl("^MT-", rownames(sce_HiPSC))

# Find ribosomal genes
is.ribosomal_HiPSC <- grepl("^RPS|^RPL",
                            rownames(sce_HiPSC))

# Calculate QC matrix for each cell,
# stored in pData os the SCEset
sce_HiPSC <- calculateQCMetrics(sce_HiPSC,
                                feature_controls = list(Rb = is.ribosomal_HiPSC,
                                                        Mt = is.mito_HiPSC))
# Plots total genes and library sizes
png("PlotQC_Mito_Ribo.png", w = 2000,
    h = 2000, res = 400)
par(mfrow = c(2, 2), cex = 1.2)
hist(sce_HiPSC$total_counts/1e+06, xlab = "Library sizes (millions)",
     main = "", breaks = 20, col = "grey80",
     ylab = "Number of cells")
hist(sce_HiPSC$total_features_by_counts, xlab = "Number of expressed genes",
     main = "", breaks = 20, col = "grey80",
     ylab = "Number of cells")
## above not workinig 
# plots reads mapped to Mt genes
hist(sce_HiPSC$pct_counts_Mt,
     xlab = "Mitochondrial proportion (%)",
     ylab = "Number of cells", breaks = 20,
     main = "", col = "grey80")

# plots reads mapped to Rb genes
hist(sce_HiPSC$pct_counts_Rb,
     xlab = "Ribosomal proportion (%)",
     ylab = "Number of cells", breaks = 20,
     main = "", col = "grey80")
dev.off()
png("AverageCount.png", w = 2000, h = 2000,
    res = 400)
# examine expression of log-means
# across all genes
ave.counts <- rowMeans(counts(sce_HiPSC))
hist(log10(ave.counts), breaks = 100,
     main = "", col = "grey80", xlab = expression(Log[10] ~
                                                    "average count"))

# Plot number of top genes
fontsize <- theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 16))


plotHighestExprs(sce_HiPSC,
       n = 30) + fontsize

dev.off()

# Plot number of cells

numcells <- nexprs(sce_HiPSC, byrow = TRUE)
png("AverageCount_SmoothScatter.png",
    w = 2000, h = 2000, res = 400)
smoothScatter(log10(ave.counts), numcells,
              xlab = expression(Log[10] ~ "average count"),
              ylab = "Number of expressing cells")
dev.off()

#add normalized log counts assay 
sce_HiPSC <- normalize(sce_HiPSC)

# PCA plot to check for potential cell
# outliers according to PC1 and PC2
# (based on general cell data)
png("PCA_cellOutliers.png", w = 2000,
    h = 2000, res = 400)
plotPCA(sce_HiPSC, by_exprs_values="logcounts") +
  fontsize
dev.off()


############## Additional data cleaning
############## steps####################################

# remove cells with low expression or
# low number of genes (lower than 3
# median absolute deviation of
# log(library size))

libsize.drop_HiPSC <- isOutlier(sce_HiPSC$total_counts,
                                nmads = 3, type = "lower", log = TRUE)
feature.drop_HiPSC <- isOutlier(sce_HiPSC$total_features_by_counts,
                                nmads = 3, type = "lower", log = TRUE)
# remove cells with high percent of
# reads mapped to Mt genes (possibly
# dead cells) (higher than 3 median
# absolute deviation)
mito.drop_HiPSC <- isOutlier(sce_HiPSC$pct_counts_Mt,
                             nmads = 3, type = "higher")
ribo.drop_HiPSC <- isOutlier(sce_HiPSC$pct_counts_Rb,
                             nmads = 3, type = "higher")
sce_HiPSC <- sce_HiPSC[, !(libsize.drop_HiPSC |
                             feature.drop_HiPSC | mito.drop_HiPSC |
                             ribo.drop_HiPSC)]

# Remove cells by Mt gene further - ISSUE 
mito.drop0.2 <- sce_HiPSC$pct_counts_feature_controls_Mt <=
  20
mito.remove0.2 <- sce_HiPSC$pct_counts_feature_controls_Mt >
  20
sce_HiPSC <- sce_HiPSC[, mito.drop0.2]

# Remove cells by Rb gene further
ribo.drop0.5 <- sce_HiPSC$pct_counts_feature_controls_Rb <=
  50
ribo.remove0.5 <- sce_HiPSC$pct_counts_feature_controls_Rb >
  50
sce_HiPSC <- sce_HiPSC[, ribo.drop0.5]

# check number of cells that genes express --didn't work
png("Cells_expresing_genes.png", w = 2000,
    h = 2000, res = 400)
numcells <- nexprs(sce_HiPSC, byrow = TRUE)
meow <- as.numeric(numcells)
hist(log2(meow), xlab = "Log2 number of cells expressing the gene",
     ylab = "Number of genes", main = "Number of cells a gene was detected")
dev.off()
# remove genes expressed in fewer than
# 1% of total cells
numcells <- nexprs(sce_HiPSC, byrow = TRUE)
genes.keep <- numcells >= 22
genes.remove <- numcells < 22
sce_HiPSC_ftGenes <- sce_HiPSC[genes.keep,
                               ]

# Plot number genes after data
# filtering
png("TopGenes_gotMapped_PostDataCleaning.png")
plotHighestExprs(sce_HiPSC,
       n = 30) + fontsize
dev.off()
# Plot number genes after data
# filtering
png("TopGenes_gotMapped_PostDataCleaning.png")
plotQC(sce_HiPSC_ftGenes, type = "highest-expression",
       n = 30) + fontsize
dev.off()

# check how many cells are removed
datRemove <- data.frame(ByLibSize = sum(libsize.drop_HiPSC),
                        ByFeature = sum(feature.drop_HiPSC),
                        GeneRemovedByCell = sum(genes.remove),
                        #ByMito1 = sum(mito.drop_HiPSC), ByMito0.2 = sum(mito.remove0.2),
                        #byRibo1 = sum(ribo.drop_HiPSC), ByRibo0.5 = sum(ribo.remove0.5),
                        CellRemaining = ncol(sce_HiPSC_ftGenes),
                        GeneRemaining = sum(genes.keep))
#note Mt and Rb genes were NOT at
# this stage

write.table(datRemove, "number_cells_genes_removed.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")
save(sce_HiPSC_ftGenes, file = "HiPSC_ftGenes_ReadyFor_ComputeSumFactors_BOC.Obj")
########## Remove Mt genes and Rb genes before
########## clustering########################
mito_ftGenes <- grep("^MT-", rownames(sce_HiPSC_ftGenes))

ribosomal_ftGenes <- grep("^RPL|^RPS",
                          rownames(sce_HiPSC_ftGenes))

mito_ribo <- c(mito_ftGenes, ribosomal_ftGenes)

sce_HiPSC_ftGenes_rmMtRb <- sce_HiPSC_ftGenes[-mito_ribo,
                                              ]
######### Normalization by deconvolution
######### method##################################
library(scran)
library(limSolve)
library(ascend)
package.version("ascend")
# Caution: takes a long time on
# laptop, better run it in cluster

sce_HiPSC_ftGenes_cp <- computeSumFactors(sce_HiPSC_ftGenes_rmMtRb,
                                          sizes = c(40, 60, 80, 100), positive = T)

# Remove zero size factors by
# converting them to the minimum size
# factor value: NORT WORKING!!!!!!!!!!!!!!!!!!!!!!!!!
Zero_sizefactor <- which(sce_HiPSC_ftGenes_cp@phenoData@data$size_factor ==
                           0)
min_size <- min(sce_HiPSC_ftGenes_cp@phenoData@data$size_factor[-Zero_sizefactor])
sce_HiPSC_ftGenes_cp@phenoData@data$size_factor[Zero_sizefactor] <- min_size

# plot the normalised data
png("sizeFactors_normalized.png")
plot(sizeFactors(sce_HiPSC_ftGenes_cp),
     sce_HiPSC_ftGenes_cp$total_counts/1e+06,
     log = "xy", ylab = "Library size (millions)",
     xlab = "Size factor")
dev.off()


# save the object before normalisation
save(sce_HiPSC_ftGenes_cp, file = "HiPSC_ftGenes_ComputeSizeFactor_before_dcvl.Obj")

# perform the KEY NORMALIZATION STEPH
sce_HiPSC_ftGenes_dcvl <- normalize(sce_HiPSC_ftGenes_cp)
# save the object after normalisation
save(sce_HiPSC_ftGenes_dcvl, file = "HiPSC_ftGenes_dcvl.Obj")


library(Seurat)
library(ascend)
library(Rtsne)
library(openxlsx)


EMSet <- EMSet(sce_HiPSC_ftGenes_cp)

EMSet <- normaliseByRLE(EMSet)
EMSet <- runPCA(EMSet,
                ngenes = 1500, 
                scaling = TRUE)

EMSet <- runCORE(EMSet,
                 conservative = FALSE,
                 remove.outliers = TRUE,
                 nres = 40,
                 dims = 10)
EMSet@clusterAnalysis
clus.info <- (EMSet@clusterAnalysis$clusters)
counts.assay <- (EMSet@assays@data$counts)
counts.assay <- as.data.frame(counts.assay)
colnames(counts.assay) <- clus.info
head(counts.assay)
row.names(counts.assay) <- make.names(row.names(counts.assay), unique=TRUE)
write.csv(counts.assay, "Day15_S1_S2.csv")
write.table(counts.assay, "d15.s1.s2.txt", sep = "\t", row.names = TRUE)
any(duplicated(rownames(counts.assay)))
cluster1_vs_all <- runDiffExpression(EMSet, 
                                     group = "cluster",
                                     condition.a = 1,
                                     condition.b = 2)
cluster2_vs_all <- runDiffExpression(EMSet, 
                                     group = "cluster",
                                     condition.a = 2,
                                     condition.b = 1)


# Create a blank workbook
OUT <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(OUT, "CLUS1")
addWorksheet(OUT, "CLUS2")

# Write the data to the sheets
writeData(OUT, sheet = "CLUS1", x = cluster1_vs_all)
writeData(OUT, sheet = "CLUS2", x = cluster2_vs_all)

# Reorder worksheets
worksheetOrder(OUT) <- c(1,2)

# Export the file
saveWorkbook(OUT, "DE_D15_Rep.xlsx")                                     


##########
one.clus <- which(EMSet@clusterAnalysis$clusters == "1")
length(one.clus)
one <- EMSet@assays@data$logcounts[,one.clus]
myh <- one["MYH7",]
myh <- as.data.frame(myh)
head(myh)
mmm <- which(myh$myh > 1)
mm <- myh[mmm,]
length(mm)

two.clus <- which(EMSet@clusterAnalysis$clusters == "2")
length(two.clus)
two <- EMSet@assays@data$logcounts[,two.clus]
myh2 <- two["MYH7",]
myh2 <- as.data.frame(myh2)
head(myh)
mmm2 <- which(myh2$myh2 > 1)
mm2 <- myh2[mmm2,]
length(mm2)

colnames(myh)


filter(one, )
subset(one, "MYH7" > 1)
EMSet <- runTSNE(EMSet, dims = 2, PCA = TRUE, seed = 1)
tsne_plot <- plotTSNE(EMSet, group = "cluster")
SeuratObject <- convert(EMSet, to = "seurat")
SeuratObject <- FindNeighbors(SeuratObject, reduction = "PCA", dims = 1:20)
SeuratObject <- FindClusters(SeuratObject, reduction="pca",dims= 1:20,k.param=30,resolution=0.01)
SeuratObject <- FindClusters(SeuratObject)
markers <- FindAllMarkers(SeuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
two <- as.data.frame(markers)
write.csv(two,"trial.markers.csv")


DimPlot(object = SeuratObject, reduction = "TSNE")
FeaturePlot(object = SeuratObject, features =c("THY1"))
## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_realData <- Rtsne(, perplexity=10, check_duplicates = FALSE)


