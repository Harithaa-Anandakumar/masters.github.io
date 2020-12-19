harimadefun <- function (height, x = "Count", color = "p.adjust", showCategory = 8, 
          font.size = 12, title = "", ...) 
{
  object <- height
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  }
  else if (x == "count" || x == "Count") {
    x <- "Count"
  }
  df <- fortify(object, showCategory = showCategory, by = x, 
                ...)
  if (colorBy %in% colnames(df)) {
    p <- ggplot(df, aes_string(x = "Description", y = x, 
                               fill = colorBy)) + theme_dose(font.size) + scale_fill_continuous(low = "red", 
                                                                                                high = "blue", name = color, guide = guide_colorbar(reverse = TRUE))
  }
  else {
    p <- ggplot(df, aes_string(x = "Description", y = x, 
                               fill = "Description")) + theme_dose(font.size) + 
      theme(legend.position = "none")
  }
  p + geom_bar(stat = "identity") + coord_flip() + ggtitle(title) +
    scale_y_reverse() +
    xlab(NULL) + ylab(NULL)+
  theme(
    panel.border = element_rect(linetype = "solid", fill = NA, size = 2, colour = "grey10"),
    panel.background = element_rect(fill = "whitesmoke"),
    panel.grid.major.y= element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = "grey30", 
                               size = 17,  hjust=1),
    axis.title.y = element_text(face = "bold", colour = "grey10",
                                size = 16)
  )
}



tiff("02_01GOenrich-1-MDC.tiff",units="in", width=10,height=14, res=300, compression = 'lzw')
harimadefun(go_enrich)
dev.off()
sto <- as.data.frame(sto)
4,1,3,2
categ <- c("two","two","two","two",rep("three",4),rep("one",8),rep("four",3))
sto$categ <- categ
sto$categ <- as.factor(sto$categ)
x <- sto$Count
a <- ifelse(sto$categ == "one", "darkred",
            ifelse(sto$categ == "two", "dimgray",
           ifelse(sto$categ == "three", "seagreen" ,"steelblue")))
      
sto$Description <- as.factor(sto$Description)
sto$Description <- factor(sto$Description, levels = unique(sto$Description))
sto$Description <- factor(sto$Description, levels = sto$Description)
p <- ggplot(sto, aes_string(x = "Description", y = x, 
                           fill = "p.adjust")) + 
  theme_dose(12) + scale_fill_continuous(low = "red", 
                                                high = "blue", 
                                                name = "p.adjust", 
                                                guide = guide_colorbar(reverse = TRUE))



p <- p + geom_bar(stat = "identity") + coord_flip()  +
  scale_y_reverse() +
  xlab(NULL) + ylab(NULL)+
  theme(
    panel.border = element_rect(linetype = "solid", fill = NA, size = 2, colour = "grey10"),
    panel.background = element_rect(fill = "whitesmoke"),
    panel.grid.major.y= element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = a, 
                               size = 20,  hjust=1),
    axis.title.y = element_text(face = "bold", colour = "grey10",
                                size = 16),
    legend.title = element_text( size = 18),
    legend.text = element_text(size = 15),
    legend.position = "left"
  )
p + scale_colour_continuous(guide = guide_legend(direction="horizontal",
                                                 label.position = "top",
                                                 label.theme = element_text(angle = 90)))
dev.off()

  
  
