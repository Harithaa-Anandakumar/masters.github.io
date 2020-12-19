#devtools::install_github("xvrdm/ggrough")
Sys.setenv("plotly_username"="Road_Runner")
Sys.setenv("plotly_api_key"="••••••••••")
library(dplyr)
library(readxl)
library(readr)
library(ggplot2)
library(reshape2)
library(viridis)
library(ggsci)
library(ggthemes)
library(ggforce)
library(purrr)
library(tidyverse)
library(DESeq2)
library(dplyr)
library(EDASeq)
library(readxl)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ggsci)
library(heatmaply)
library(ggrough)
dat <- read_csv("barPlotInfo_STAT_all_CT_5.txt")
dat <- as.data.frame(dat)
pwd
dat <- as.data.frame(read.table("data/barPlotInfo_STAT_all_CT_5.txt",sep="\t",
                                header=TRUE))
dat <- as.data.frame(read.table("barPlotInfo_all_CT_5.txt",sep="\t",
                                header=TRUE))
dat <- as.data.frame(read.table("data/decon/barPlotInfo_ge_CT_5.txt",sep="\t",
                                header=TRUE, stringsAsFactors = FALSE))

tnr <- as.data.frame(read.table("barPlotInfo_STAT_sp_CT_5.txt",sep="\t",
                                header=TRUE))
dim(dat)
row.names(dat) <- dat$X
dat <- dat[,-1]
head(dat)
dat.prop <- dat/colSums(dat)

dat.prop[1:3,1:3]
maxab <- apply(dat, 1, max)

head(maxab)

n1 <- names(which(maxab > 5.00))

dat1 <- dat[-which(row.names(dat) %in% n1), ]
colnames(dat) <- c("X","Sample1","Sample2","Sample3")
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)

heatmap(as.matrix(dat),Rowv = NA, Colv = NA, col = scaleyellowred)

library(pheatmap)
pheatmap(as.matrix(dat),
         scale = "none",
         #breaks= max(8),
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="BrBG")))(100))
dev.off()
colnames(dat) <- c("X", sample)
melted <- melt(dat)
melted$X <- factor(melted$X,levels = unique(melted$X))
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")
x <- ggplot(melted,aes(x = variable, y = X)) + 
  geom_point(aes(size = value, fill = X), alpha = 0.75, shape = 21)+  
  scale_size_continuous(limits = c(0.000000, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        #legend.title = element_blank(),
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")  
x + guides(fill=FALSE)  
ggsave("bacterial_abundance.png")
#scale_fill_manual(values = colours, guide = FALSE) + 
  #scale_y_discrete(limits = rev(levels(melted$variable)))

p <- ggplot(melted,aes(x = variable, y = X)) + 
  geom_point(aes(size = value, fill = X), alpha = 0.75, shape = 21)

ggrough(p,options)


dim(dat1)
dat1$x <- rownames(dat1)
melted <- melt(dat1)
xrange <- range(melted$y)
melted$y <- rep(1:27)
yrange <- range(melted$value)
xkcdaxis(xrange,yrange)
p <-ggplot(melted, (aes(x=x, y=value, fill=variable)))+
  geom_bar(stat="identity") +
  theme_minimal()+
  scale_fill_tableau()

ggplot(dat1,(aes(x=)))
p+scale_fill_manual(values= mypal(19))

mypal <- colorRampPalette(brewer.pal(6, "Set1"))

theme(panel.background = element_rect(fill="#223583"))


options <- list(
  Background=list(roughness=8),
  GeomCol=list(fill_style="dots",  angle_noise=0.5, fill_weight=2))
get_rough_chart(p, options)
get_rough_chart(p)

dat2 <- dat2[-1,]
dat2 <- as.data.frame(t(dat1))
colnames(dat2) <- dat$X
rownames(dat2) 
sample <- c("D1","D2","D3")  
total <- c(37300000,30800000,34600000)
human <- as.numeric(c(24600000,20300000,21600000))  
bacteria <-  as.numeric(c(0.001308,0.000942,0.000532))
rna <- as.numeric(c(1687691,1704332,3222209))
viral <- c(0.0002978,0.0002188,0.0002570)
h_mil <- human/10^6
bac <- bacteria * 10 ^6
##number of bacterial reads per million human mapped reads
brpm <- bac/h_mil
##number of viral reads per million human mapped reads
vrpm <- (viral*10^6)/(human/10^6)
vrpm
### rRNA 
rpm <- data.frame(brpm,vrpm)
rpm$sample <- sample
colnames(rpm) <- c("Bacteria","Virus","Sample")
mrpm <- melt(rpm)

ggplot(mrpm, aes(x=Sample, y=value, colour=variable))+
  geom_point(size=5)+
  labs(y="Viral / Bacterial reads per million human mapped reads")+
  theme_bw()+
  theme(text=element_text(size=15))
ggsave("contaminants_permill_human_reads.png")

bar <- data.frame(sample,human,bacteria,viral)
rownames(bar) <- bar$sample
bar <- bar[,-1]
bar.prop <- (bar/rowSums(bar))*100
as.integer(bar.prop)
bar.prop <- as.data.frame(bar.prop) 
bar.prop$sample <- rownames(bar)
bar.melted <- melt(bar.prop)

#!!!!!worked!!!!this-factor-numeric-conv-imp!!!!!
dat2[] <- lapply(dat2, function(x) as.numeric(as.character(x)))


p <- ggplot(bar.melted, (aes(x=sample, y=value, fill=variable)))+
  geom_bar(stat="identity") +
  theme_minimal()+
  scale_fill_tableau()
p +scale_y_log10()
v <- (human/total)*100
h_percent <- (human/total)*100
v_percent <- as.numeric(viral/total)*100
#v_percent <- format(v_percent,scientific = F)
as.numeric(v_percent)
b_percent <- (bacteria/total)*100
m <- h_percent + v_percent + b_percent + rRNA_percent
rRNA_percent <- c(4.52,5.54,9.31)
other_per <- 100 - m
a=3.506702e-09/100
b=7.4634656/100
c <- b_percent/100
class(d)
d <-as.numeric(as.character(dat2$Acinetobacter))
percent_of_percent <- dat2*c           
rownames(percent_of_percent) <- sample
bar <- data.frame(h_percent,v_percent,other_per)
bar
bar$sample <- sample
bar$rna <- rRNA_percent
new_bar <- bar
new_bar$bacteria <- b_percent
###DOT PLOT with only 4 groups - human, rrna, viral and bacterial
breaks1 = c(60,30,5.0,0,1.264185e-9)
melt_new <- melt(new_bar)
melt_new$group <- c(rep("Human",3),rep("Viral",3),rep("Multi-Mapped",3),rep("rRNA",3),rep("Bacterial",3))
melt_new$group <- factor(melt_new$group, levels = c("Human","Multi-Mapped","rRNA","Bacterial","Viral"))
ggplot(melt_new, (aes(x=sample, y=(value),colour =group)))+
  geom_point(stat="identity",na.rm=TRUE,size=3) +
  scale_color_nejm()+
  labs(y="% contribution")+
  theme_bw()+
  labs(face="bold")+
  scale_y_log10(breaks = breaks1, labels = function(x) format(x,scientific=TRUE))+
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 9),
        axis.title.x = element_text(colour = "black", face = "bold", size = 11),
        axis.title.y = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")
ggsave("combine_dotplot.png")  

percent_of_percent$sample <- sample
merged <- merge(bar,percent_of_percent, by="sample")
mel.bar <- melt(merged)
mel.bar <- melt(bar)
mel.bar$log <- log10(mel.bar$value)
##Go through each row and determine if a value is zero
row_sub = apply(mel.bar, 1, function(row) all(row != "0.000000e+00" ))
##Subset as usual
mel.bar <- mel.bar[row_sub,]
dim(mel.bar)
mel.bar$variable 
newbar <- mel.bar[mel.bar$new=="bacteria",]\
class(newbar)
newbar <- newbar[-(1:6),]
mel.bar$group <- c(rep("human",3),rep("viral",3),rep("multi-mapped",3),rep("rRNA",3),rep("bacteria",72))
library(plotly)
library(scales)

#this weirdly works - don't ask how
ggplot(mel.bar, (aes(x=sample, y=value,fill =variable)))+
  geom_bar(position= "fill",stat="identity",na.rm=TRUE) +
  scale_fill_viridis_d()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw()
##lets try new shit here
ggplot(mel.bar, (aes(x=sample, y=value,fill =variable)))+
  geom_bar(stat="identity",na.rm=TRUE,position="dodge") +
  coord_trans(y="log10")+
  scale_fill_viridis_d()

#re-order factor levels in mel.bar$group to change their legend labelling
mel.bar$group <- factor(mel.bar$group, levels = c("human","multi-mapped","rRNA","viral","bacteria"))
##this dot chart weirdly works -- gotta add ticks though

breaks = 0.000001**(0.000000001:100 * 0.01)
breaks = c(60,30,5.0,0,3.240405e-11,1.264185e-10)
ggplot(mel.bar, (aes(x=sample, y=(value),colour =group)))+
  geom_point(stat="identity",na.rm=TRUE) +
  scale_color_lancet()+
  labs(y="% contribution")+
  #scale_y_continuous(trans="log10",labels =scales::comma)
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
               # labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(breaks = breaks, labels = function(x) format(x,scientific=TRUE))+
  theme_bw()
ggsave("dotplot.png")  
#coord_trans(y="log10") +
  #scale_y_continuous(trans="log10",labels=function(x) format(x,scientific=FALSE))
#annotation_logticks(scaled=FALSE)

scale_y_continuous()


ggplot(mel.bar,aes(x=sample,y=value))+
  geom_point(aes(colour=variable))

p + scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))


min <- min(mel.bar$value)
max <- max(mel.bar$value)
ggplotly(p)


library(reshape2)
dcast(dat, X)
genus <- dat$X
x <- as.data.frame(t(dat))
colnames(x) <- genus
xdim(x)
x <- x[-1,]
x
rownames(x) <- bar$sample
x <- x[,-1]
x
bar
class(x)
7.983914e-10 * 7.427746e-10
x <- as.numeric(x)
x <- as.matrix(x)
7.6722338/b_percent * 100
b_percent/x[,]
b_percent/x
class(b_percent)

b_percent/x[,-1]

1000/x[,-1]

x <- as_tibble(x)
x %>%
  mutate_all(funs(changed = 2.00/.))
x <- as.data.frame(x)
head(x)
x$new <- c(3,3,3)
colnames(x)


plswork <- x %>% 
  mutate_at(vars(-new), funs(. / new))
plswork
a <- 3.506702e-09 /dat[,2]
b <- 3.058442e-09/dat[,3]
c <- 1.537572e-09/dat[,4]
new <- data.frame(a,b,c)
new$genus <- dat$X
colnames(new) <- c ("D1","D2","D3","X")
tail(plswork)

meh <- as.data.frame(t(new))
colnames(meh) <- dat$X
rownames(meh)
meh <- meh[-4,]
meh

rownames(bar) <- bar$sample
bar <- bar[,-5]
merged <- merge(meh,bar,by=0)
row.names(merged) <- merged$Row.names
merged <- merged[,-1]
melt.merg <- melt(merged,id=Row.names) 

dontknow <- as.data.frame(t(merged))
dontknow <- dontknow[-1,]
colnames(dontknow) <- c ("D1","D2","D3")
row.names(dontknow)
dontknow
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
merged <- as_tibble(merged)
merged$Row.names <- as.factor(merged$Row.names)
merged$h_percent <- as.numeric.factor(merged$h_percent)
merged$v_percent <- as.numeric.factor(merged$v_percent)
worked <- melt(z, id = "Row.names")
meh
na.omit(meh)
meh <-
ggplot(worked, (aes(x=Row.names, y=value, fill=variable)))+
  geom_bar(stat="identity") +
  theme_minimal()+
  scale_fill_tableau()
dev.off()
p + scale_y_log10()

y<- as.data.frame(t(merged))

z <- as.data.frame(t(y))

z
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)
data

worked %>% 
  filter_all(all_vars(!is.infinite(.)))

worked$value <- wait
 <- worked$value
wait <- as.integer(wait)
mm <- worked[apply(worked, 1, function(x) all(is.finite(x))), ]
mm             
