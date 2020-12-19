library(tidyverse)
library(broom)
library(clusterProfiler)

dds <- DESeqDataSetFromMatrix(countData = raw_counts_df, colData = coldat, design=~ new)
dds <- dds[ rowSums(counts(dds)) > 5, ]
combined_thesis <- assay(dds)
colnames(combined_thesis) <- dds$meh
combined_thesis <- as.data.frame(combined_thesis)
combined_thesis$Gene <- rownames(combined_thesis)
combined_thesis<- cbind(combined_thesis$Gene, combined_thesis)
rownames(combined_thesis) <- NULL
combined_thesis <- combined_thesis[,-53]
colnames(combined_thesis) <- c("Gene", dds$meh)
write.csv(combined_thesis,"combined_forCIB.csv")
write.table(combined_thesis, "combined_forCIBERSORT_rawCounts.txt", sep="\t", quote = F)
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
ntop=5000
Pvars <- rowVars(assay(onlyCM))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(Pvars)))]
PCA <- prcomp(t(assay(onlyCM)[select, ]), scale = F)
plot(PCA$x[,1],PCA$x[,2])
vsd <- vst(dds)
pc <- plotPCA(onlyCM, intgroup="new", returnData=TRUE)
pc
onlyCM <- dds[,dds$Group =="CM"]
vsd <- vst(dds)
onlyCM <- vst(onlyCM)
#checking how much variation princ component 1 accounts for to get a sense of how meaningful these clusters are;
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
colour_list <- list(
  new = c(Adult ="#80B1D3", Fetal ="#FB8072", Rh = "#BEBADA",
            CM = "#FFFFB3", EHM = "#8DD3C7"))
tiff("02_03PCA_MDC.tiff",units="in", width=7,height=5, res=300, compression = 'lzw') 
ggplot(data=pc,aes(x=PC1,y=PC2, color=new)) +
  geom_point(size=5)+
 # scale_color_manual(values =c("Adult" ="darkred", "Fetal" ="orangered", 
  #                             "Rh" = "#BEBADA",
  #                             "CM" = "darkolivegreen", "EHM" = "salmon") )+
  #scale_color_tableau(palette="Color Blind")+
  scale_color_few()+
  xlab(paste0("PC1: ", percentVar[1],"% variance"))+
  ylab(paste0("PC2: ", percentVar[2], "% variance"))+
  theme(
    legend.title = element_blank(),
    legend.background = element_rect(fill="whitesmoke"),
    legend.direction = "horizontal",
    legend.text = element_text( size=14, 
                                face="bold"),
    legend.position = c(0.6, 0.9),
  
    panel.border = element_rect(linetype = "solid", fill = NA, size = 2, colour = "grey10"),
    panel.background = element_rect(fill = "whitesmoke"),
    panel.grid.major= element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(face = "bold", color = "grey30", 
                               size = 12, hjust=1),
    axis.text.y = element_text(face = "bold", color = "grey30", 
                               size = 12,  hjust=1),
    axis.title.y = element_text(face = "bold", colour = "grey10",
                                size = 16),
    axis.title.x = element_text(face = "bold", colour = "grey10",
                                size = 16)
  )
dev.off()
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
  dplyr::slice(1:500) %>%
  # extract the column (gene name) from the table
  pull(column) %>%
  # retain unique gene names only
  unique()

gene_loadings <- PCA$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::filter(gene %in% top_genes) 


loadings_plot <- ggplot(gene_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_label_repel(data = gene_loadings, 
                   aes(x = PC1, y = PC2, label = gene))



plotPCA(onlyCM, intgroup="new")
new <- assay(onlyCM)[rownames(assay(onlyCM)) %in% three,]
colnames(new) <- vsd$new
new <- new[apply(new, MARGIN = 1, FUN = function(x) sd(x) != 0),]
tiff("99_01heatmap-glycolysis_onlyCM.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
an_col <- as.data.frame(onlyCM$new)
rownames(an_col) <- onlyCM$x
colnames(an_col) <- c("group")


colour_list <- list(
                    group = c(Adult ="#80B1D3", Fetal ="#FB8072", Rh = "#BEBADA",
                              CM = "#FFFFB3", EHM = "#8DD3C7"),
                    cluster = c("1" ="#F4CAE4","2"= "#CBD5E8" ,
                                "3"="#FDCDAC","4"= "#B3E2CD" ))
new <- assay(onlyCM)[rownames(assay(onlyCM)) %in% top_genes,]
colnames(new) <- onlyCM$x
colorRampPalette(rev(brewer.pal(n = 5, name = "Accent")))(8)

clusDesignation <- as.data.frame(cutree(as.hclust(out$tree_row), 4))
colnames(clusDesignation) <- c("cluster")
clusDesignation$genes <- rownames(clusDesignation)
clus <- as.data.frame(clusDesignation[order(clusDesignation$cluster),])
write.csv(clus, "trialMarkers1.csv")
summary(out)

tiff("02_02heatmap-MDC.tiff",units="in", width=6,height=9, res=300, compression = 'lzw') 
pheatmap(na.omit(new),
         scale="row",
         #color=magma(20),
         #main = "",
         #cutree_cols = 7,
         cluster_rows=TRUE,
         clustering_distance_rows="correlation",
         cex=1,
         clustering_distance_cols="manhattan", 
         clustering_method="complete", 
         border_color=FALSE,
         cluster_cols = TRUE,
         #color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(length(breaksList)),
         #color=meh,
         color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdGy")))(1000),
         annotation_col= an_col,
         annotation_row = clusDesignation,
         annotation_colors = colour_list,
         cutree_rows = 4,
         cutree_cols = 3,
         show_rownames = FALSE,
         show_colnames = FALSE,
       # breaks = mat_breaks,
         treeheight_row = 0)
dev.off()
# calulate the correlations
?cor
colnames(new) <- onlyCM$meh
r <- cor(new, method="pearson")
round(r,2)
# Compute a matrix of correlation p-values
p.mat <- cor_pmat(new)
library(ggcorrplot)
ggcorrplot(r)

########THIS IS GOOD!!!##########
ggcorrplot(r, 
           hc.order = TRUE, 
           #type = "lower",
           p.mat = p.mat,
           lab = FALSE)

#########################

ggcorrplot(s)
s <- r[rownames(r)=="Adult_Heart",]
t <- as.data.frame(colMeans(s))
t$group <- colnames(s)
t$group <- as.factor(t$group)
colnames(t) <-c("corr", "group")

u <- t[t$group!= "Adult_Heart",]
u <- u[u$group!= "CM-Rh",]

p <- ggplot(u, aes(group, corr))
p+geom_boxplot()

colnames(new) <- onlyCM$meh



# Using median
tiff("02_01PearsonCorr-MDC.tiff",units="in", width=5,height=5, res=300, compression = 'lzw') 
p <- u %>%
  mutate(class = fct_reorder(group, corr, .fun='median')) %>%
  ggplot( aes(x=reorder(group, corr), y=corr, fill=group)) + 
  geom_boxplot(size=1) +
  ylab("Pearson Corelation\n") +
  theme(legend.position="none") +
  scale_fill_tableau(palette =  "Superfishel Stone") +
  xlab("")

p + theme(
  panel.border = element_rect(linetype = "solid", fill = NA, size = 2, colour = "grey10"),
  panel.background = element_rect(fill = "whitesmoke"),
  panel.grid.major.y= element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(face = "bold", color = "grey30", 
                             size = 12, angle = 45, hjust=1),
  axis.text.y = element_text(face = "bold", color = "grey30", 
                             size = 12,  hjust=1),
  axis.title.y = element_text(face = "bold", colour = "grey10",
                              size = 16)
  )
dev.off()
##########################################################
normalized_counts <- counts(dds, normalized=TRUE)

dds <- DESeq(dds)
plotDispEsts(dds)
res1_unshrunken <- results(dds, contrast = c('new','CM','Adult'), alpha = 0.05, lfcThreshold = 0.58)
res_tableOE <- lfcShrink(dds,contrast = c('new','CM','Adult'), res=res1_unshrunken)

summary(res1_unshrunken)
### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
res_tableOE_tb <- res2_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sigOE <- res_tableOE_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

## Order results by padj values
top20_sigOE_genes <- res_tableOE_tb %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20) 		#Extract the first 20 genes
top2000_sigOE_genes <- res_tableOE_tb %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=5000) 
### write it out onto a file
CM_vs_EHM <- res_tableOE_tb %>% 
  arrange(padj)
write.csv(CM_vs_EHM, "CM_vs_EHM.csv")
# Create tibbles including row names
mov10_meta <- coldat %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

## normalized counts for top 20 significant genes
top20_sigOE_norm <- normalized_counts %>%
  filter(gene %in% top20_sigOE_genes)
# Gathering the columns to have normalized counts to a single column
gathered_top20_sigOE <- top20_sigOE_norm %>%
  gather(colnames(top20_sigOE_norm)[2:52], key = "samplename", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top20_sigOE)


gathered_top20_sigOE <- inner_join(mov10_meta, gathered_top20_sigOE)

## plot using ggplot2
ggplot(gathered_top20_sigOE) +
  geom_point(aes(x = gene, y = normalized_counts, color = new)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))
### Extract normalized expression for significant genes from the OE and control samples (4:9), and set the gene column (1) to row names
norm_OEsig <- normalized_counts[,c(1,2:52)] %>% 
  filter(gene %in% top2000_sigOE_genes) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 
### Annotate our heatmap (optional)
annotation <- mov10_meta %>% 
  select(meh, new, samplename) %>% 
  data.frame(row.names = "samplename")

annotation <- mov10_meta %>% 
  select(meh, new)
annotation <- as.data.frame(annotation)
rownames(annotation) <- make.names(dds$new, unique=T)

### Set a color palette
heat_colors <- brewer.pal(9, "YlOrRd")


colnames(norm_OEsig) <- dds$new
### Run pheatmap

pheatmap(norm_OEsig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)


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




#######################
# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="MYOM1", intgroup="new", returnData=TRUE)

# Plotting the MOV10 normalized counts, using the samplenames (rownames of d as labels)
ggplot(d, aes(x = new, y = count, color = new)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
 # geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  #ggtitle("MYL2") +
  theme(plot.title = element_text(hjust = 0.5))


#############################
dat <- new[1:52,]
dat
dat <- normalized_counts %>% 
        filter(gene %in% three[1:50])

gathered <- dat %>%
  gather(colnames(dat)[2:52], key = "samplename", value = "normalized_counts")

joined <- inner_join(mov10_meta, gathered)
joined <- joined %>%
  filter(new != "Fib")
joined <- joined %>% filter( new != "Rh")
xb <- joined %>%
  group_by(new, gene) %>%
  dplyr::summarize(median = median(normalized_counts, na.rm=TRUE),
                   sd = sd(normalized_counts, na.rm=TRUE),
                   mean = mean(normalized_counts, na.rm=TRUE))



ggplot() +
  geom_point(data=joined, aes(x = gene, y = normalized_counts, color=new),  alpha=0.3) +
  geom_point(data=xb,aes(x=gene, y=median, colour=new, size=0.5),  alpha=1)+
  #scale_shape_manual(values=c(21,22,24,25))+
 geom_errorbar(data=xb, mapping = aes(x = gene, y = median,  ymin = mean - sd, ymax = mean + sd, color=new),size=0.2, width=.2)+
  xlab("Genes") +
  scale_y_log10() +
  #ylab("log10 Normalized Counts") +
  #ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  scale_colour_tableau(palette = "Tableau 10")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15))

g +  scale_color_tableau(guide = guide_legend(keywidth = 5,
                                             keyheight = 5))+
  guides(colour=guide_legend(override.aes = list(size = 4)))

dev.off()

library(scales)
show_col(tableau_color_pal(palette="Tableau 10")(9))
#########################################################
library(org.Hs.eg.db)
require(clusterProfiler)
data(geneList)
gene <- clusDesignation$genes
univer <- rownames(onlyCM)
gene.df <- bitr(univer, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

ok <- merge(gene.df, clusDesignation, by.x="SYMBOL",by.y="genes")
ok = ok[order(ok$cluster),]


num.val <- PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "1" |  PC == "2") %>%
  # for each PC
  group_by(PC) %>%
  # sort descending value
  dplyr::arrange(desc(abs(value))) %>%
  dplyr::filter(column %in% top_genes)

d <- merge(ok, num.val, by.x="SYMBOL", by.y="column")

d = read.csv(your_csv_file)
## assume 1st column is ID
## 2nd column is FC
d = read_xlsx("data/DE_D15_Rep.xlsx")

## feature 1: numeric vector
geneList = d[,6]
## feature 2: named vector
names(geneList) = as.character(d[,4])
## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)



library(enrichplot)
hari <- ok$ENTREZID
names(hari) <- ok$cluster
four <- ok %>%
  dplyr::filter(cluster == "4")
myList <- list(one$ENTREZID,two$ENTREZID,three$ENTREZID,four$ENTREZID)
names(myList) <- c("one",'two',"three","four")
go_enrich <- enrichGO(gene = two$ENTREZID,
                      universe= gene.df$ENTREZID, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "CC",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

upsetplot(go_enrich)
upsetplot(ck)
ck@compareClusterResult$Description

go_enrich@result$Description 

go_enrich <- filter(go_enrich, Description %in% monkey)



barplot(go_enrich, 
        drop = TRUE, 
        showCategory =50,
        title = "GO Biological Pathways",
        font.size = 14)

dotplot(go_enrich,
        font.size=14)
dev.off()
class(monkey)

library(clusterProfiler.dplyr)

monkey <- c("contractile fiber","contractile fiber part","sarcomere","I band","A band","M band",
             "actin filament bundle","Z disc","myosin complex")
emapplot(go_enrich)


intwp <- c("L-type voltage-gated calcium channel complex","I band","contractile actin filament bundle",
  "neuron to neuron synapse","neurotransmitter receptor complex","presynapse","contractile fiber","contractile fiber part","sarcomere","I band","A band","M band",
  "actin filament bundle","Z disc","myosin complex","presynaptic membrane","integral component of synaptic membrane",
  "presynaptic active zone","postsynaptic membrane","dendritic spine","GABA-ergic synapse","synaptic cleft")


intwp <- c("synaptic cleft","presynapse")
goplot(go_enrich, showCategory = 10)
dotplot(go_enrich, showCategory=50)

ck <- compareCluster(geneClusters = myList,
                     OrgDb = org.Hs.eg.db,
                     universe = gene.df$ENTREZID,
                     fun = "enrichGO")
head(as.data.frame(ck))
dotplot(ck)
dev.off()
(summary(ck))
gsecc <- gseGO(geneList=geneList,
               ont="CC",
               pvalueCutoff = 0.9,
               OrgDb=org.Hs.eg.db, 
               verbose=F)


require(DOSE)
x = enrichDO(one$ENTREZID, qvalueCutoff=1, pvalueCutoff=1)
## convert gene ID to Symbol
edox <- setReadable(go_enrich, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'

p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
p1

heatplot(edox, foldChange = geneList)
go_enrich <- enrichGO(gene = three$ENTREZID,
                      universe= gene.df$ENTREZID, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "CC",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
secondClus <- go_enrich@result[1:10,]
firstClus <- go_enrich@result
thridClus <- go_enrich@result[1:4,]
fourthClus <- go_enrich@result[1:4,]
4,1,3,2
firstClus <- firstClus[firstClus$Description!="myosin complex",]
fourthClus <- fourthClus[fourthClus$Description!="collagen trimer",]
fourthClus <- fourthClus[fourthClus$Description!="extracellllar matrix component",]
secondClus <- secondClus[5:10,]
secondClus <- secondClus[2:5,]
secondClus <- secondClus[secondClus$Description!="collagen-containing extracellular matrix",]
sto <- rbind(secondClus,thridClus,firstClus,fourthClus)
sto
