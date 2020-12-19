
#Figure 4.9
write.csv(gen, "data/genesFromAtas.csv")
gen <- read.csv("data/genesFromAtas.csv")
# perform pca using prcomp func, choosing the top n number of genes 
ntop=2000
Pvars <- rowVars(assay(vsd))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]

PCA <- prcomp(t(assay(vsd)[select, ]), scale = T)


top_genesHM <- PCA %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  dplyr::filter(PC == "1" |  PC == "2") %>%
  # for each PC
  group_by(PC) %>%
  # sort descending value
  dplyr::arrange(desc(abs(value))) %>%
  # take top 5 rows of each PC group
  dplyr::slice(1:2000) %>%
  # extract the column (gene name) from the table
  pull(column) %>%
  # retain unique gene names only
  unique()


mew <- as.data.frame(assay(vsd)) %>%
  as_tibble(rownames = "gene") %>%
  filter(gene %in% top_genesHM)

mew <-  as.data.frame(mew)
rownames(mew) <- mew$gene
mew <- mew[,-1]
colnames(mew) <- paste(vsd$group, vsd$specification, sep = " - " )
# Pairwise correlation between samples (columns)
cols.cor <- cor(mew, 
                use = "pairwise.complete.obs", 
                method = "pearson")

# Pairwise correlation between rows (genes)
rows.cor <- cor(t(mew), 
                use = "pairwise.complete.obs", 
                method = "pearson")

#rownames(mew) <- NULL
# Plot the heatmap


# Pairwise correlation between samples (columns)
forhc <- mew
colnames(forhc) <- paste(vsd$group, vsd$specification, sep = " - " )
cols.cor1 <- cor((forhc), 
                 use = "pairwise.complete.obs", 
                 method = "pearson")

library(dendsort)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
hc.complete = hclust(dist(t(mew)), method="complete")
hc.complete$label <- paste(vsd$group, vsd$specification, sep = " - " )
hc.complete.rows = hclust(dist((mew)), method="complete")


find_coordinates = function(n, gaps, m = 1:n){
  if(length(gaps) == 0){
    return(list(coord = unit(m / n, "npc"), size = unit(1 / n, "npc") ))
  }
  
  if(max(gaps) > n){
    stop("Gaps do not match with matrix size")
  }
  
  size = (1 / n) * (unit(1, "npc") - length(gaps) * unit("4", "bigpts"))
  
  gaps2 = apply(sapply(gaps, function(gap, x){x > gap}, m), 1, sum) 
  coord = m * size + (gaps2 * unit("4", "bigpts"))
  
  return(list(coord = coord, size = size))
}

draw_colnames = function(coln, gaps, vjust_col, hjust_col, angle_col, ...){
  coord = find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = vjust_col, hjust = hjust_col, rot = angle_col, gp = gpar(...))
  
  return(res)
}

assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames",
  ns = asNamespace("pheatmap")
)

tiff("heatMapThesis.tiff",units="in", width=10,height=10, res=300, compression = 'lzw') 
pheatmap(
  mew,
  scale = "row",
  cluster_cols=hc.complete,
  cluster_rows =hc.complete.rows,
  clustering_distance_cols = as.dist(1 - cols.cor),
  clustering_distance_rows = as.dist(1 - rows.cor),
  cutree_cols = 6,
  cutree_rows = 6,
  annotation_colors = heatmapColScale,
  #annotation_col = colAn,
  #annotation_row = rowAn,
  #show_rownames = FALSE,
  #show_colnames = FALSE,
  labels_col = hc.complete$label,
  treeheight_row = 0,
  angle_col = 315,
  fontsize_col = 7,
  col=inferno(100)
  #  col = inferno(length(breaksList)),
  # breaks = breaksList
)

newForCorr <- newForCorr %>%
  as_tibble(rownames = "gene") %>%
  select("gene","CM","EHM","Fetal_Heart","Adult_Heart")
library(ggplotify)

x <- as.ggplot(pheatmap(newForCorr[2:5],
         show_rownames = FALSE,
         cluster_cols = FALSE,
         scale="row",
         color=viridis(100)))
x | huh
newForCorr <- sapply(split(seq_len(ncol(forCorr)),colnames(forCorr)),function(cis) rowMeans(forCorr[,cis,drop=F]))

corrplot(cols.cor1, 
         method = "color",
         tl.col='grey30',
         addrect = 3,
         order = 'hclust',
         tl.cex = 0.5,
         number.cex = .7,
         col=inferno(200),
         #type="upper",
         #  bg = "black",
         diag = FALSE) 

forCorr <- mew
colnames(forCorr) <- vsd$group
r <- cor(forCorr, method="pearson")
round(r,2)

s <- r[rownames(r)=="Adult_Heart",]
t <- as.data.frame(colMeans(s))
t$group <- colnames(s)
t$group <- as.factor(t$group)
colnames(t) <-c("corr", "group")

u <- t[t$group!= "Adult_Heart",]

# Using median
tiff("PearsonCorr-thesis.tiff",units="in", width=4,height=4, res=300, compression = 'lzw') 
huh <- u %>%
  mutate(class = fct_reorder(group, corr, .fun='median')) %>%
  ggplot( aes(x=reorder(group, corr), y=corr)) + 
  geom_boxplot(aes(fill = stage(group, after_scale = alpha(fill, 0.9))))+
  # geom_boxplot(fill=after_stat(),alpha=0.4) +
  ylab("Pearson Corelation\n") +
  theme(legend.position="none") +
  scale_fill_manual(values = colscale)+
  #scale_fill_tableau(palette =  "Superfishel Stone") +
  xlab("")+
  theme_few()+
  guides(fill=FALSE)




dev.off()

library(ggplotify)
################################################
onnu <- as.data.frame(assay(vsd)) %>%
  as_tibble(rownames = "gene") %>%
  filter(gene %in% top_genesHM)
onnu <-  as.data.frame(onnu)
rownames(onnu) <- onnu$gene
onnu <- onnu[,-1]
colnames(onnu) <- vsd$group
onnuPh <- sapply(split(seq_len(ncol(onnu)),colnames(onnu)),function(cis) rowMeans(onnu[,cis,drop=F]))
onnuPh <- onnuPh %>%
  as_tibble(rownames = "gene") %>%
  select("gene","CM","EHM","Fetal_Heart","Adult_Heart")
ph1000 <- as.ggplot(pheatmap(onnuPh[2:5],
                        show_rownames = FALSE,
                        cluster_cols = FALSE,
                        treeheight_row = 0,
                        #scale="row",
                        color=viridis(100),
                        fontsize = 15,
                        angle_col = 45))


colnames(onnu) <- vsd$group
r <- cor(onnu, method="pearson")
round(r,2)

s <- r[rownames(r)=="Adult_Heart",]
t <- as.data.frame(colMeans(s))
t$group <- colnames(s)
t$group <- as.factor(t$group)
colnames(t) <-c("corr", "group")

u <- t[t$group!= "Adult_Heart",]

cor1500 <- u %>%
  mutate(class = fct_reorder(group, corr, .fun='median')) %>%
  ggplot( aes(x=reorder(group, corr), y=corr)) + 
  geom_boxplot(aes(fill = stage(group, after_scale = alpha(fill, 0.9))))+
  # geom_boxplot(fill=after_stat(),alpha=0.4) +
  ylab("Pearson Corelation\n") +
  theme(legend.position="none") +
  scale_fill_manual(values = colscale)+
  #scale_fill_tableau(palette =  "Superfishel Stone") +
  xlab("")+
  theme_few()+
  scale_y_continuous(limits=c(0.20,0.80), breaks = seq(0.20,0.80,.10))+
  theme(axis.text.y = element_text(size=15, colour="black"),
        axis.text.x = element_text(color="black", angle=45, hjust = 1, face="bold",size=15))+
  guides(fill=FALSE)


rendu <- as.data.frame(assay(vsd)) %>%
  as_tibble(rownames = "gene") %>%
  filter(gene %in% gen$x)
rendu <-  as.data.frame(rendu)

rownames(rendu) <- rendu$gene
rendu <- rendu[,-1]
colnames(rendu) <- vsd$group
renduPh <- sapply(split(seq_len(ncol(rendu)),colnames(rendu)),function(cis) rowMeans(rendu[,cis,drop=F]))
renduPh <- renduPh %>%
  as_tibble(rownames = "gene") %>%
  select("gene","CM","EHM","Fetal_Heart","Adult_Heart")
ph380 <- as.ggplot(pheatmap(renduPh[2:5],
                             show_rownames = FALSE,
                             cluster_cols = FALSE,
                            treeheight_row = 0,
                             #scale="row",
                             color=viridis(100),
                            fontsize = 15,
                            angle_col = 45))


colnames(rendu) <- vsd$group
r <- cor(rendu, method="pearson")
round(r,2)

s <- r[rownames(r)=="Adult_Heart",]
t <- as.data.frame(colMeans(s))
t$group <- colnames(s)
t$group <- as.factor(t$group)
colnames(t) <-c("corr", "group")

u <- t[t$group!= "Adult_Heart",]

cor380 <- u %>%
  mutate(class = fct_reorder(group, corr, .fun='median')) %>%
  ggplot( aes(x=reorder(group, corr), y=corr)) + 
  geom_boxplot(aes(fill = stage(group, after_scale = alpha(fill, 0.9))))+
  # geom_boxplot(fill=after_stat(),alpha=0.4) +
  ylab("Pearson Corelation\n") +
  theme(legend.position="none") +
  scale_fill_manual(values = colscale)+
  #scale_fill_tableau(palette =  "Superfishel Stone") +
  xlab("")+
  theme_few()+
  scale_y_continuous(limits=c(0.20,0.80), breaks = seq(0.20,0.80,.10))+
  theme(axis.text.y = element_text(size=15, colour = "black"),
        axis.text.x = element_text(color="black", angle=45, hjust = 1, face="bold", size=15))+
  guides(fill=FALSE)

#patchwork magick!

p1 <- (ph1000 / cor1500 )+
  plot_layout(heights = c(2,1))

p2 <- (ph380 / cor380) +
  plot_layout(heights = c(2,1))

tiff("01Corr.tiff",units="in", width=6,height=6, res=300, compression = 'lzw') 
cor1500
dev.off()

library(patchwork)
patchwork <- (p1 |  p2) 

patchwork[[1]] <- patchwork[[1]] + plot_layout(tag_level = 'new')
patchwork[[2]] <- patchwork[[2]] + plot_layout(tag_level = 'new')

tiff("multiCorr-thesis.tiff",units="in", width=6,height=6, res=300, compression = 'lzw') 
patchwork + plot_annotation(tag_levels = c('A', '1'),
                            tag_sep = '.')  
dev.off()


############################

eig <- get_eigenvalue(PCA)
eig <- eig[1:min(8, nrow(eig)), , drop = FALSE]
trying <- seq(1:8)
text_labels <- round(eig, 1)
df.eig <- as.data.frame(eig)
df.eig$dimension <- seq(1:8) 

ggplot(df.eig, aes(dimension, cumulative.variance.percent))+
  geom_bar(stat = "identity", fill="beige", colour="grey")+ 
  geom_line(color = "grey") + 
  geom_point(shape = 19, color = "red")+
  geom_text(label = text_labels$cumulative.variance.percent, vjust = -0.8, 
                       hjust = 1, size=4.5)+
  scale_x_continuous(breaks=seq(1,8,1))+
  labs(y="Percentage of Cummulative Variance Explained",
       x="Principal Components")+
  theme_few()+
  theme(axis.text.y = element_text(size=12, colour = "black"),
        axis.text.x = element_text(color="black"))


####################################

