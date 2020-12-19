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
library(gridExtra)
library(grid)
library(PCAtools)
f_files<- list.files("/Users/bhuvaneswari/Desktop/Thesis/malte/skeletal_muscle", pattern = "*.txt", full.names = T)


read_in_feature_counts<- function(file){
  cnt<- read_tsv(file, col_names =F, comment = "#")
  #cnt<- cnt %>% dplyr::select(-Chr, -Start, -End, -Strand, -Length)
  return(cnt)
}
raw_counts<- map(f_files, read_in_feature_counts)
raw_counts_df<- purrr::reduce(raw_counts, inner_join,by = "X1") 
colnames(raw_counts_df) <- c("Geneid","SRR4296453",
                             "SRR4296454",
                             "SRR4296455",
                             "SRR4296456",
                             "SRR4296457",
                             "SRR4296458",
                             "SRR4296459",
                             "SRR4296460",
                             "SRR4296461",
                             "SRR4296462")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- raw_counts_df$Geneid
G_list <- getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id","hgnc_symbol","external_gene_name"),
  values=genes,
  mart= mart,
  uniqueRows = TRUE)

G_list <- data.frame(
  genes[match(G_list$ensembl_gene_id,genes)],
  G_list)
colnames(G_list) <- c(
  "genes",
  c("ensembl_gene_id", "gene_biotype", "external_gene_name"))

G_list
raw_counts_df <- merge(raw_counts_df,G_list,by.x="Geneid",by.y="ensembl_gene_id")
#anno<- left_join(raw_counts_df,G_list,by=c(raw_counts_df$Geneid == G_list$ensemble_gene_id))
row.names(raw_counts_df) <- make.names(raw_counts_df$external_gene_name, unique = TRUE)
raw_counts_df <- raw_counts_df[,-1]
colnames(raw_counts_df)
raw_counts_df <- raw_counts_df[,-(11:13)]

write.table(raw_counts_df, "hicks_nature_combined_counts.txt", sep ="\t", quote = FALSE)

tomerge <- read.table("hicks_nature_combined_counts.txt")
f_files<- list.files("/Users/bhuvaneswari/Desktop/Thesis/malte/skeletal_muscle", pattern = "*txt", full.names = T)


read_in_feature_counts<- function(file){
  cnt<- read_tsv(file, col_names =T, comment = "#")
  cnt<- cnt %>% dplyr::select(-Chr, -Start, -End, -Strand)
  return(cnt)
}
raw_counts<- map(f_files, read_in_feature_counts)
raw_counts_df<- purrr::reduce(raw_counts, inner_join) 
raw_counts_df <- as.data.frame(raw_counts_df)
rownames(raw_counts_df) <- raw_counts_df$Geneid
raw_counts_df <- raw_counts_df[,-1]
raw_counts_df <- merge(raw_counts_df,tomerge,by=0)
cnt <- raw_counts_df
cnt <- cnt[,-(1:2)]

rownames(length) <- length$gene
raw_counts_df <- merge(cnt,length,by=0)
row.names(raw_counts_df) <- raw_counts_df$Row.names
colnames(raw_counts_df)
new <- raw_counts_df[,34:35]
cnt <- raw_counts_df[,-1]
cnt <- cnt[,-(33:34)]
##calculating TPM
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

okay <- cnt$p786sESMcontrol7w3_xx_totalrna_sr_Tiburcy_J_1.bam / new$raw_counts_df.Length
p786sESMcontrol7w3_xx_totalrna_sr_Tiburcy_J_1.bam <- okay/sum(okay, na.rm = TRUE) * 1e6

tpm(cnt$p786sESMcontrol7w3_xx_totalrna_sr_Tiburcy_J_1.bam, new$raw_counts_df.Length)
tpms <- apply(cnt,2,function(x) tpm(x,new$raw_counts_df.Length))
rownames(tpms) <- rownames(cnt)
tpms <- as.data.frame(tpms)
tpms$p786sESMcontrol7w3_xx_totalrna_sr_Tiburcy_J_1.bam <- p786sESMcontrol7w3_xx_totalrna_sr_Tiburcy_J_1.bam


write.csv(tpms, "TPM-normalized-modified2_skeletal_muscle.csv")
tpms <- tpms[-1,]
class(tpms)
tpms <- rbind(spec, tpms)
tpms <- rbind(group, tpms)
group <- as.vector(coldat$group)
spec <- as.vector(coldat$batch)
rownames(raw_counts_df) <- raw_counts_df$Row.names
raw_counts_df <- raw_counts_df[,-1]
raw_counts_df <- na.omit(raw_counts_df)

write.table(raw_counts_df, "skeletal_muscle_proj_rawCounts.txt", sep="\t",quote=FALSE)
coldat <- colnames(raw_counts_df)
write.csv(coldat,"coldat_skeletal_muscle.csv")

coldat <- read.csv("/Users/bhuvaneswari/Desktop/Thesis/malte/skeletal_muscle/coldat_skeletal_muscle.csv")

rownames(coldat) <- coldat[,1]
coldat$team <- as.factor(coldat$team)
coldat$batch <- as.factor(coldat$batch)
coldat$group <- as.factor(coldat$group)

cnt <- read.table("/Users/bhuvaneswari/Desktop/Thesis/malte/skeletal_muscle/rawCounts_skeletal_muscle.txt")


##########
dds <- DESeqDataSetFromMatrix(countData = raw_counts_df, colData = coldat, design=~ batch )
dds <- dds[ rowSums(counts(dds)) > 5, ]
vsd <- vst(dds)
before_batch_correction <- plotPCA(vsd, intgroup="team")
assay(vsd) <- limma::removeBatchEffect(assay(vsd), dds$team)
after_batch_correction <- plotPCA(vsd, intgroup="team")


one <- plotPCA(vsd, intgroup="team")
p <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)

tiff("pca_skeletal_muscle",units="in", width=10,height=5, res=300, compression = 'lzw') 
biplot(p,
       lab=NULL,
       colby = 'group',
       pointSize = 5,
       legendPosition = 'left', legendLabSize = 13, legendIconSize = 6.0,
       shape = 'team', shapekey = c('malte'=15, 'hicks'=17))
dev.off()

ntop=5000
Pvars <- rowVars(assay(vsd))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]                   
PCA <- prcomp(t(assay(vsd)))

PCA <- prcomp(t(assay(vsd)[select, ]), scale = F)

meh <- na.omit(assay(vsd)[select,])
colnames(meh) <- vsd$group
tiff("heatmap_skeletal_muscle",units="in", width=10,height=7, res=300, compression = 'lzw')
pheatmap(meh,
         treeheight_row = 0,
         scale="row",
         main = "heatmap-top 5000 variant genes",
         cutree_cols = 3,
         cluster_rows=TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="BrBG")))(100),
         show_rownames = FALSE,
         show_colnames = TRUE)
dev.off()
#install.packages("corrplot")
library(corrplot)
m <- cor(meh)
library(RColorBrewer)

tiff("corr-plot_skeletal_muscle",units="in", width=10,height=10, res=300, compression = 'lzw') 
corrplot(m, order="FPC", col = brewer.pal(n = 8, name = "RdBu"),
         tl.col='black')
dev.off()


