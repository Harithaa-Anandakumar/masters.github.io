library(EnhancedVolcano)
library(dplyr)
##read the Orthology table 
orth.tab <- read.table("rhesus/mart_export.txt", sep ="\t")

##retain useful columns
orth.tab <- orth.tab %>%
  dplyr::filter(V7 == "1") %>%
  distinct(V1,V6)

##read in human files 
f_files<- list.files("data/subgroup", pattern = "*.txt", full.names = T)
f_files
read_in_feature_counts<- function(file){
  cnt<- read_tsv(file, col_names = T, comment = "#")
  cnt<- cnt %>% dplyr::select(-Chr,-Start, -End, -Strand)
  return(cnt)
}
raw_counts<- map(f_files, read_in_feature_counts)
raw_counts_df<- purrr::reduce(raw_counts, inner_join) 
raw_counts_df <- data.frame(raw_counts_df)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- raw_counts_df$Geneid
######
G_list <- getBM(
  filters= "external_gene_name", 
  attributes= c("ensembl_gene_id","hgnc_symbol","external_gene_name"),
  values=genes,
  mart= mart,
  uniqueRows = TRUE)

rownames(G_list) <- make.names(G_list$ensembl_gene_id, unique=TRUE)
rownames(orth.tab) <- make.names(orth.tab$V1, unique = TRUE)
## first get the ensembleIDs for the human gene names
raw_counts_df <- merge(raw_counts_df,G_list,by.x="Geneid",by.y="external_gene_name")

##now merge the ensembleIDs from both the orth.table and our counts table
one <- merge(orth.tab, G_list, by=0)
one.one <- merge(one, raw_counts_df, by.x ="V1", by.y="ensembl_gene_id")

##convert human gene names to the 1:1 orthologous gene names of the rhesus
rownames(one.one) <- make.names(one.one$V6,unique=TRUE)
colnames(one.one)

##Remove columns that are not needed 
one.one <- one.one[,-(1:7)]
one.one <- one.one[,-12]
#have a 00_identifier in front of all human samples
colnames(one.one) <- paste("00",colnames(one.one),sep="_")

## read in monkey files
f_files<- list.files("rhesus/counts", pattern = "*.txt", full.names = T)
f_files
read_in_feature_counts<- function(file){
  cnt<- read_tsv(file, col_names = T, comment = "#")
  cnt<- cnt %>% dplyr::select(-Chr,-Start, -End, -Strand)
  return(cnt)
}
raw_counts<- map(f_files, read_in_feature_counts)
raw_counts_df<- purrr::reduce(raw_counts, inner_join) 
raw_counts_df <- data.frame(raw_counts_df)
raw_counts_df$Geneid <- gsub("gene-","\\1",raw_counts_df$Geneid)
rownames(raw_counts_df) <- raw_counts_df$Geneid
raw_counts_df
raw_counts_df1 <- raw_counts_df1[,-1]
##merge the monkey samples with the human samples based on the common monkey gene names
one.three <- merge(one.one, raw_counts_df, by=0)

rownames(one.three) <- one.three$Row.names
colnames(one.three)
##we're going to create a matrix of the length sizes of each gene for each sample -- which can be used to normalize
#same genes in the different species are of diff lengths 
no.length <- one.three$`00_Length`
monkey.length <- one.three$Length

##now again clean up data
one.three <- one.three[,-(1:2)]
one.three <- one.three[,-(11:12)]

#make length matrix
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

norm.mat <- rep.col(no.length, 10)
norm.mat <- cbind(norm.mat,monkey.length,monkey.length
                  )
colnames(norm.mat) <- NULL
##write the length matrix out into a file
write.table(norm.mat, "lengthMat_RH.txt", sep="\t", quote=FALSE)
write.table(one.three, "rawCounts_RH.txt", sep="\t", quote=FALSE)

coldat <- read.csv("data/colDat_RH.csv")
rownames(coldat) <- coldat$x
coldat$Specification <- as.factor(coldat$Specification)
coldat$Folder <- as.factor(coldat$Folder)
coldat$Group <- as.factor(coldat$Group)
coldat$type <- as.factor(coldat$type)
one <- read.table("rhesus/rawCounts_RH.txt", sep="\t")
two <- read.table("rhesus/lengthMat_RH.txt", sep="\t")
dds <- DESeqDataSetFromMatrix(countData = one, colData = coldat, design=~Group)


three <- read_xlsx("Desktop/Thesis/self_contained_thesis/mmc5_Tib.xlsx")
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
colnames(one)
human <- one[,1:10]
monk <- one[,11:12]
hum_length <- two$V1
monk_length <- two$V12
tpm_monk <- apply(monk,2,function(x) tpm(x,monk_length))
tpm_hum <- apply(human,2,function(x) tpm(x, hum_length))
colnames(tpm_hum)
tpmNorm_RH <- merge(tpm_hum, tpm_monk, by= 0)
write.csv(tpmNorm_RH,"Desktop/Thesis/self_contained_thesis/rhesus/tpmNorm_RH.csv")
## add matrix into this slot in a DESEQ2 object -- it'll now use this as normalization factor
assays(dds)[["avgTxLength"]] <- two
vsd <- vst(dds)
assay(vsd)

plotPCA(vsd, intgroup="Group")

dev.off()
# perform pca using prcomp func, choosing the top n number of genes 
ntop=2000
Pvars <- rowVars(assay(vsd))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]

PCA <- prcomp(t(assay(vsd)[select, ]), scale = T)



l <- PCA$x
#l <- plotPCA(vsd, intgroup=c("group","specification"), returnData=TRUE)
#percent variation each component accounts for
percentVarp <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
tiff("rhesusPCA.tiff",units="in", width=5,height=3, res=300, compression = 'lzw') 
ggplot(as.data.frame(l), aes(PC1, PC2, color = vsd$Group))+
  geom_point(size =4) +
  #ggtitle("PCA-using all genes") +
  theme(plot.title = element_text(hjust=0, size=8))+
  #geom_label_repel(data = as.data.frame(l),aes(PC1, PC4),
  #             label = vsd$specification) +
  xlab(paste0("PC1: ", round(percentVarp[1],digits=2),"% variance"))+
  ylab(paste0("PC2: ", round(percentVarp[2],digits = 2), "% variance"))+
  guides(shape=FALSE)+
  scale_color_tableau(palette = "Tableau 10")+
#  scale_color_manual(values = colscale)+
  theme_few()+
  theme(legend.title = element_blank())

dev.off()
norm_counts_RH <- assay(vsd)
write.csv(norm_counts_RH,"normCounts_RH.csv")
dds <- DESeq(dds)
plotCounts(dds, "TNNI1", intgroup = "Group")
res1 <- results(dds, contrast = c("Group","Cardiomyocytes","RH-CM"))
head(result.df)
# Remove genes with no p-value
uninformative_genes <- is.na(result.df$padj) 
result.df <- result.df[!uninformative_genes,]

res1 <- lfcShrink(dds,
          contrast = c('Group','Cardiomyocytes','RH-CM'), res=res1, type ="normal")
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
# perform pca using prcomp func, choosing the top n number of genes 
ntop=2000
Pvars <- rowVars(assay(vsd))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]

PCA <- prcomp(t(assay(vsd)[select, ]), scale = T)
l <- PCA$x
#l <- plotPCA(vsd, intgroup=c("group","specification"), returnData=TRUE)
#percent variation each component accounts for
percentVarp <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
ggplot(as.data.frame(l), aes(PC1, PC3, color = vsd$Group))+
  geom_point(size =3) +
  #ggtitle("PCA-using all genes") +
  theme(plot.title = element_text(hjust=0, size=8))+
  #geom_label_repel(data = as.data.frame(l),aes(PC1, PC4),
  #             label = vsd$specification) +
  #xlab(paste0("PC1: ", round(percentVarp[1],digits=2),"% variance"))+
  ylab(paste0("PC2: ", round(percentVarp[3],digits = 2), "% variance"))+
  guides(shape=FALSE)+
 # scale_color_manual(values = colscale)+
  theme_few()+
  theme(legend.title = element_blank())
top_genes3 <- PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "3"| PC == "1") %>%
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
library(factoextra)
filtered <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  #dplyr::filter(gene %in% top_genes3) %>%
  dplyr::filter(PC3 < -0.01  )  %>%
  #group_by(PC) %>%
  # sort descending value
  dplyr::arrange(desc(abs(PC3)))
  #dplyr::slice(1:100) 
write_csv(filtered, "RH_PCA_drivingDown.csv")
fviz_pca_biplot( PCA,
                 axes = c(1, 3),
                 geom.ind = "point",
                 fill.ind = vsd$Group,
                 col.ind =  "lightgrey",
                 pointsize = 2, pointshape=21,
                 #palette = colscale,
                 addEllipses = FALSE,
                 # col.ind = vsd$group,
                 #variables
                 alpha.var = 0.1,
                 col.var = "contrib",
                 select.var = list(names= filtered$gene),
                 gradient.cols = pal,
                 #c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,
                 labelsize = 3, 
                 axes.linetype=NA)  + ggtitle("") +
  labs(x= 'PC1', y= 'PC3')+
  guides(fill=FALSE)+
  theme_few()



new <- assay(vsd)[select,]

colnames(new) <- dds$Group
new <- new[apply(new, MARGIN = 1, FUN = function(x) sd(x) != 0),]


genes <- c("MYH7","MYBPC3","ACTN2","TNNI3","MYL2","ATP2A2","MYH6",
           "TNNI1","MYL7","NKX2-5","GATA4")
new <- assay(vsd)[rownames(assay(vsd)) %in% three,]
colnames(new) <- dds$Specification
meh <- vsd[, vsd$type == "CM"]
new <- assay(meh)[rownames(assay(meh)) %in% three,]
colnames(new) <- meh$Specification
tiff("heatmap_usingGeneList_RH.tiff",units="in", width=5,height=6, res=300, compression = 'lzw') 

pheatmap(na.omit(new),
         scale="row",
         color=viridis(20),
         #main = "",
         #cutree_cols = 7,
         cluster_rows=TRUE,
         cluster_cols = TRUE,
         #color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         show_rownames = FALSE,
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

atrialMarkers <- c("NPPB","HAMP","MYBPHL","NPPA","PLA2GA","COMP","TCEAL2","SLP1",
                   "DHRS9","HP","HEY1","NR2F1","NR2F2",
                   "TBX5","ATP2A2","MYH6","GNAO1","KCNA5","KCNJ1","
                        KCNK3","KCNQ2","CACNA1D","HOXA3","PDGFRA","PDE8B")

ventricularMarkers <- c("DLK1","IRX4","MYL2","XDH","TMEM190","HYAL2",
                        "CPNE4","CYP1A1","IRX5","C3orf23",
                        "HAND1","HEY2","MYH7","GJA1","NAV1","KCNJ2","KCNJ4",
                        "LPL","COL12A1","SCUBE3")

cell_surfaceMarkers <- c("VCAM","TMEM71","TMEM173","TMEM151A","ROR2","CDH2",
                         "GPR177","NCAM1","CSCR7","TMEM66","GPC2","LGR2","
                         PDGFRA","FZD4","CD99",
                         "SIRPA","GPR37","FLRT2","EPCAM","ANPEP","LEPREL1")

mesoderm <- c("EOMES","FOXF1","GATA6","GSC","MESP1","MIXL1","WNT5A")

early_cardiac_commitment <- c("GATA4","HOXB2","IRX3","GATA5","MEIS2",
                              "IRX5","MEIS1","HEY1","ZNF503","ZBTB16",
                              "PBX3","RARB","RXRA","ZFPM1","PPARG",
                              "RUNX1T1")
early_cardiacProgenitor_EnrichedTFs <- c("GATA4","KLF11","SMAD3","HAND2",
                                         "LHX2","NR2F1","TSHZ2","ZFHX3",
                                         "ZFPM2")
progenitor_cm_regulatory <- c("SMAD6","SOX18","ISL1","MEF2C","TCF21","NKX2-5",
                              "TBX18","KLF2","ZBTB20","MEF2D","WT1","SMARDCD3","TBX2")

progenitor_cm_functional <- c("MYL2","KCNA5","KCNH2","MYL7","CAMK2B","CACNB2",
                              "DES","ATP2A2","MYL9","TCAP","MYL4","CACNA1C","MYH6","TNNI1",
                              "TNNT2","RYR2","TPM2","TTN")



fa <- c("ACADM","CPT2","CPT1B","ETFDH","HADHA","HADHB","MAPK14")

oxre <- c("HIBADH","IMPDH2","NDUFA10","NDUFA11","NDUFA12","NDUFB6","NDUFS1","NDUFV1",
          "ACAD9","ACADM","ADHFE1","ALDH4A1","ALDH5A1",
          "ALDH6A1","ALDH7A1","AASS","AJFM1","CAT","CRYZ","CYP27A1","CYBSA","CYB5A","COX681",
          "DHTKD1","DHFR","DLD","ETFA","ETFB","ETFDH",
          "FTH1","GLUD1","GCDH","GPD2","GRHPR","HADH","HSD17B10","IDH2","MDH2","OGDH","PRDX3",
          "PCYOX1","PYROXD2","PDHA1","PDHB","SLC25A12","SLC25A13","SUOX","TMLHE","UQCRF51","UQCRFS1")

glycolysis <- c("BPGM","TPI1","ALDOA","ALDOC","ENO1","ENO2","GPI","GAPDH","HK2","LDHA","PFKL",
                "PFKP","PGK1","PGAM1")

up_human_heart <- c("CASQ2","MB","MYOM2","TCAP","MYH11","TNNI3","S100A1","DES",
                    "HRC","MYOM1")

up_hh_lateCMs <- c("MYH7","MYL2","TNNI3K","HSPB7","PLN","CSRP3","ACTN2","RBM20",
                   "TRIM63","CORIN")

up_earlyCMs <- c("NKX2-5","IRX4","TBX2","COL2A1",
                 "ISL1","HAND1","ID2","LEF1","IRS1",
                 "MDK")
fiftyGenes <- c('Myh6', 'Atp2a2', 'Pln', 'Cox6a2', 'Cox7b', 'Ndufa1', 'Uqcrq', 'Atp5e', 'Cox7a1', 'Cox6c', 'Tnni3', 'Fabp3', 'Myl2', 'Pgam1', 'Tubb5', 'Nme1', 'Gm5506', 'Eif5a', 'Ngfrap1', 'Cks1b', 'Cdkn1c', 'Mest', 'Gpc3', 'H2afz', 'Tnni1', 'Mif', 'Hmgn2', 'Gyg', 'Myl7', 'Myl4', 'Hadha', 'Ryr2', 'Srl', 'Nfib', 'Nfia', 'Klf6', 'Itm2b', 'Ech1', 'Phyh', 'Oxct1', 'Gpc1', 'Fhl2', 'Mt1', 'Mgst3', 'Acadl', 'Lpl', 'Brp44l', 'D830015G02Rik', 'Lars2', 'S100a1')

sarcomeric_genes <- c("MYBPC3","MYH7","MYL2","MYL3","ACTC1","TNNC1",
                      "TNNI3","TNNT2","TPM1","ACTN2","CSRP3","MYOZ2",
                      "TCAP","TTN",
                      "CAPN6","TNNT3","SCIN","MYO15B","MYO3B","FLIP1L",
                      "CAPN14","TLN2","SGCG",
                      "SLCA4A3","MYO7B","SGCA","SGCD","SSPN","NRAP","TRIM63","MYO19B","CAPN3",
                      "ITGB1BP2","MYH7B","KLHL41","TMOD1","FILIP1","DMD","MYPN")


fiftyGenes <- toupper(fiftyGenes)
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


