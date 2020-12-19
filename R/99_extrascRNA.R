my_packages <- c("tidyverse", "broom", "coefplot", "cowplot","gapminder", "GGally", "ggrepel", "ggridges", "gridExtra","here", "interplot", "margins", "maps", "mapproj","mapdata", "MASS", "quantreg", "rlang", "scales", "survey","srvyr", "viridis", "viridisLite", "devtools")
install.packages(cellrangerRkit, repos = "http://cran.rstudio.com")
library(Seurat)
library(dplyr)
library(cowplot)

data <- read.table(file = "filtered_spots_count_matrix_all_weeks.tsv")
dim(data)
data[1:30,1:3]
object.size(data)
object.size(as.matrix(data))
data <- as.matrix(data)
pbmc <- CreateSeuratObject(counts = data, min.cells = 3, min.features  = 200, project = "10X_PBMC", assay = "RNA")
pbmc
# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
str(pbmc)
mito.genes <- grep(pattern = "^MT-", x = rownames(pbmc@assays[["RNA"]]), value = TRUE)

percent.mito <- Matrix::colSums(pbmc@assays[["RNA"]][mito.genes, ])/Matrix::colSums(pbmc@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
#Seurat v2 function, but shows compatibility in Seurat v3
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito") 
#in case the above function does not work simply do:
pbmc$percent.mito <- percent.mito

VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mito")
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito >  -Inf & percent.mito < 0.05 )


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc.new <- FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
#pbmc.new
#head(x = HVFInfo(object = pbmc.new))
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCounts_RNA", "percent.mito"))
pbmc <- RunPCA(object = pbmc,  npcs = 20, verbose = FALSE)
DimPlot(object = pbmc, reduction = "pca")
FeaturePlot(object = pbm, features = "TSPAN6")
VariableFeaturePlot(object = pbm)
dev.off()
####
length(pbmc@)
PCAPlot(pbmc, dim.1 = 1, dim.2 =2, do.label=TRUE)




######
DimHeatmap(object = pbm, reduction = "pca", cells = 200, balanced = TRUE)
ElbowPlot(object = pbm)
?ElbowPlot
ElbowPlot(new, ndims = 20, reduction = "pca")
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)
pbmc <- FindClusters(pbmc, reduction="pca",dims= 1:20,k.param=30,resolution=0)
new <- RunPCA(object = new,  npcs = 30, verbose = FALSE)
new <- FindClusters(new, reduction.type="pca",dims.use= 1:20,k.param=30,resolution=2.1)
new <- RunTSNE(new)
DimPlot(object = new, reduction = "tsne")
pbmc <- RunTSNE(object = pbmc, dims.use = 1:20, do.fast = TRUE)
DimPlot(object = pbmc, reduction = "tsne")
?RunICA
pbmc <- RunICA(pbmc,nics =7)
pbmc <- FindVariableFeatures(object = pbmc, 
                             mean.function = ExpMean, 
                             dispersion.function = LogVMR,
                             x.low.cutoff = 0.2, 
                             x.high.cutoff = 10, 
                             y.cutoff = 0.5)  # if this fails, experiment with the num.bin setting
length(pbmc@)
genes.dissoc <- c("ATF3", "BTG2", "CEBPB", "CEBPD", "CXCL3", "CXCL2", "CXCL1", "DNAJA1", "DNAJB1", "DUSP1", "EGR1", "FOS", "FOSB", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1A", "HSPA1B", "HSPA8", "HSPB1", "HSPE1", "HSPH1", "ID3", "IER2", "JUN", "JUNB", "JUND", "MT1X", "NFKBIA", "NR4A1", "PPP1R15A", "SOCS3", "ZFP36")
#### seurat <- ?(?, genes.list = list(?), ctrl.size = 20, enrich.name = "genes_dissoc")
pbmc <- AddModuleScore(pbmc, genes.list = list(genes.dissoc), ctrl.size = 20, enrich.name = "genes_dissoc")
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tail(pbmc.markers)
pbmc(assay)
hmm <- (pbmc$seurat_clusters)
dim(hmm)
write.csv(pbmc, "workpls")
write_sp


current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7,8,9)
new.cluster.ids <- c("Compact Ventricular Myocardium","Trabecular ventricular myocardium","Trabecular ventricular myocardium","Trabecular ventricular myocardium",
                     "Atrial Myocardium","Outflow tract/large vessels","Atrioventricular mesenchyme & valves","mediastinal mesenchyme & vessels",
                     "cavities with blood","epicardium")

pbmc@active.ident <- plyr::mapvalues(x = pbmc@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(object = pbmc, reduction = "tsne", do.label = TRUE, pt.size = 0.5)
pbmc <- StashIdent(object = pbmc, save.name = "ClusterNames")
pbmc$ClusterNames_1
new$seurat_clusters
data_to_write <- as.data.frame(as.matrix(new$seurat_clusters))
Cells(new)
sc <- Bis
saveRDS(new,"new.rds")
tail(data_to_write)
str(pbmc)
export_data
rna[1:3,1:3]
head(pbmc@meta.data)
pbmc[1:6,1:6]

new1 <- RenameCells(pbmc, new.names = paste0(pbmc$ClusterNames,"_", Cells(pbmc)))
one <- GetAssayData(new1, assay = "RNA", slot = "counts")
one <- as(Class="matrix", object = one)
one %>% mutate(col = str_remove(col, '_X.*'))
colnames(one) <- sub("_X.*","",colnames(one))

library(dplyr)
library(stringr)
str(new1)
head(one)
dim(one)
one <- as.data.frame(one)
three <- sub("_X*","",one[,1:1995])
two <- for ( col in 1:ncol(one)){
    colnames(one)[col] <-  sub("_X*", "", colnames(one)[col])
  }

two <- as.data.frame(pbmc.markers)
write.csv(two,"seeifitworks.csv")
one1 <- GetAssayData(new1, assay = "RNA", slot = "data")
one <- as(Class="matrix", object = one)
write.csv(one,"workpls.csv")
