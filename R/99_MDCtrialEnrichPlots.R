library(org.Hs.eg.db)
require(clusterProfiler)
library(enrichplot)
library(clusterProfiler.dplyr)
d = read_xlsx("data/DE_D15_Rep.xlsx")
d2 = read_xlsx("data/DE_D15_Rep.xlsx", sheet=2)

geneList <- as.vector(d$foldChange)
geneList2 <- as.vector(d2$foldChange)
neu <- merge(d, gene.df, by.x="id", by.y="SYMBOL")
neu2 <- merge(d2, gene.df2, by.x="id", by.y="SYMBOL")
geneList2 <- as.vector(neu2$foldChange)
names(geneList) <- as.character(gene.df$ENTREZID)
names(geneList2) <- as.character(gene.df2$ENTREZID)
geneList2 <- sort(geneList2, decreasing = TRUE)
univer <- as.data.frame(d[,1])
univer2 <- as.data.frame(d2[,1])
gene.df2 <- bitr(univer2$id, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

gene.df <- bitr(univer$id, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
balh <- rownames(raw_counts_df)
balhUniv <- bitr(balh, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

go_enrich <- enrichGO(gene = gene.df$ENTREZID,
                      universe= balhUniv$ENTREZID, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "CC",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
heatplot(go_enrich)
dev.off()
go_enrich2 <- enrichGO(gene = gene.df2$ENTREZID,
                      universe= balhUniv$ENTREZID, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "CC",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
dotplot(go_enrich, showCategory=30)
go_enrich@result$Description

neu <- merge(d, gene.df, by.x="id", by.y="SYMBOL")

geneList <- as.vector(neu$foldChange)
names(geneList) <- as.character(neu$ENTREZID)
geneList <- sort(geneList, decreasing = TRUE)


ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

ego2 <- gseGO(geneList     = geneList2,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

monkey <- c("sarcomere","system development","cardiovascular system developmetn",
            "cell differentiation","extracellular matrix",
            "extracellular matrix organization",
            "focal adhesion","cell migration","cardiac muscle tissue development",
            "muscle contraction","cardiac muscle contraction","heart development")

ego <- filter(ego3, 
                    Description %in% monkey)
ego2 <- filter(ego2, 
              Description %in% monkey)

dotplot(ego2, showCategory=30)

ego@result$Description

  
upsetplot(ego)
ggo <- groupGO(gene     = gene.df$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

ggo <- groupGO(gene     = gene.df2$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)
myList <- list(gene.df$ENTREZID,gene.df2$ENTREZID)
names(myList) <- c("one",'two')
ck <- compareCluster(geneClusters = myList,
                     OrgDb = org.Hs.eg.db,
                     universe = balhUniv$ENTREZID,
                     fun = "enrichGO")
dotplot(ck)
heatplot(ego2)
ck@compareClusterResult$Description %in% monkey
ckfil<- filter(ck, 
               Description %in% monkey)
