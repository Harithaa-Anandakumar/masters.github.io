if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("circlize")


library(ComplexHeatmap)
library(circlize)
library(cluster)


colnames(top_genes_exprs.mx) <- dds$Specification
mat <- top_genes_exprs.mx[,-(11:12)]
colnames(mat)
mat <- mat[,-17]
out <- 
rowan <- as.data.frame(heat[,1:2])
rowan <- as.data.frame(rowan[,-1])
rownames(rowan) <- heat$Row.names
pheatmap(na.omit(new), 
                show_rownames=F, cluster_cols=T, cluster_rows=T, scale="row",
                cex=1, clustering_distance_rows="euclidean", cex=1,
                clustering_distance_cols="euclidean", 
                clustering_method="complete", border_color=FALSE,
         annotation_row = rowan)

clusDesignation <- as.data.frame(cutree(as.hclust(out$tree_row), 5))
colnames(clusDesignation) <- c("cluster")
summary(out)
head(new)
heat <- merge(clusDesignation, new, by=0)




rownames(heat) <- heat$Row.names
colnames(heat)
#Set annotation
ColAnn <- data.frame(colnames(heat))
colnames(ColAnn) <- c("Sample")
ColAnn <- HeatmapAnnotation(df=ColAnn, which="col")
myCol <- colorRampPalette(c("navyblue", "white", "red"))(100)
myBreaks <- seq(-2,2, length.out=100)
RowAnn <- data.frame(heat$cluster)
colnames(RowAnn) <- c("Cluster")
colours <- list("Cluster"=c("1"="royalblue","2"="red3", "3"="yellow", "4"="green","5"="black","6"="blue","7"="pink"))
RowAnn <- HeatmapAnnotation(df=RowAnn, col=colours, which="row")
colnames(heat)
boxAnnCol <- HeatmapAnnotation(boxplot=anno_boxplot(heat[,3:27], border=FALSE, gp=gpar(fill="#CCCCCC"), lim=NULL, pch=".", size=unit(2, "mm"), axis=FALSE, 
                                                    annotation_width=unit(c(1, 7.5),"cm")))

boxAnnRow <- rowAnnotation(boxplot=row_anno_boxplot(heat[,3:27], border=FALSE, gp=gpar(fill="#CCCCCC"), lim=NULL, pch=".", size=unit(3, "cm"), axis=FALSE,
                                                    annotation_width=unit(c(3), "cm")))
                         
monkey <- as.matrix(heat[,3:27])
library(colorRamps)
hmap <- Heatmap(monkey,
                name="Transcript Z-score",
                col=colorRamp2(myBreaks, myCol),
                heatmap_legend_param=list(color_bar="continuous", legend_direction="horizontal", legend_width=unit(5,"cm"), title_position="topcenter", title_gp=gpar(fontsize=15, fontface="bold")),
                
                #Split heatmap rows by gene family
                split=heat$cluster,
                
                #Row annotation configurations
                cluster_rows=TRUE,
                show_row_dend=FALSE,
                
                #row_title="Transcript", #overridden by 'split' it seems
                row_title_side="left",
                row_title_gp=gpar(fontsize=15, fontface="bold"),
                show_row_names=FALSE,
                row_names_side="left",
                row_title_rot=0,
                
                #Column annotation configuratiions
                cluster_columns=TRUE,
                show_column_dend=TRUE,
                column_title="Samples",
                column_title_side="top",
                column_title_gp=gpar(fontsize=15, fontface="bold"),
                column_title_rot=0,
                show_column_names=TRUE,
                
                #Dendrogram configurations: columns
                clustering_distance_columns=function(x) as.dist(1-cor(t(x))),
                clustering_method_columns="ward.D2",
                column_dend_height=unit(30,"mm"),
                
                #Dendrogram configurations: rows
                clustering_distance_rows=function(x) as.dist(1-cor(t(x))),
                clustering_method_rows="ward.D2",
                row_dend_width=unit(30,"mm"))
                
                #Annotations (row annotation must be added with 'draw' function, below)
                #top_annotation_height=unit(0.5,"cm"),
                #top_annotation=ColAnn,
                
                #bottom_annotation_height=unit(3, "cm"),
                #bottom_annotation=boxAnnCol)




draw(hmap + RowAnn + boxAnnRow, heatmap_legend_side="left", annotation_legend_side="right")

clusDesignation$cluster <- as.factor(clusDesignation$cluster)
clusDesignation$gene <- rownames(clusDesignation)
trailMarkers <- clusDesignation[order(clusDesignation$cluster),]
write.csv(trailMarkers, "trialMarkers.csv")

