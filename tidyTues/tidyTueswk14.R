devtools::install_github("thebioengineer/tidytuesdayR")

brewing_materials <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-03-31/brewing_materials.csv')
beer_taxed <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-03-31/beer_taxed.csv')
brewer_size <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-03-31/brewer_size.csv')
beer_states <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-03-31/beer_states.csv')
library(lubridate)
library(tidytuesdayR)
library(ggtextures)
library(xkcd)
library(numform)
library(patternplot)
library(ggtextures)
library(magick)
library(ggimage)
library(ggplot2)
library(png)
library(grid)
library(ggimage)
library(tidyr)
library(dplyr)
library(ggthemes)
library(stringr)
library(ggsci)
library(ggrough)
library(reshape2)

class(beer_states$year)
beer_states <- as.data.frame(beer_states)
beer_states$year <- as.Date(as.character(beer_states$year),format="%Y")

beer_states$state <- as.factor(beer_states$state)

huh <- beer_states %>% 
  group_by(year)%>%
  dplyr::filter(state == "total")

beer_states %>% 
  group_by(year)%>%
  dplyr::filter(state == "total") %>%
  ggplot(aes(x=year, y=barrels, color=type)) +
  geom_line(size=3)+
  scale_x_date(date_labels = "%Y")+
  scale_y_continuous(label = ff_denom()) +
  theme_pander() +
  scale_color_rickandmorty()

##this DOES not work
brewing_materials$date <- paste(brewing_materials$year, brewing_materials$month, sep="-") %>% 
  ymd() %>% 
  as.Date()
#####
brewing_materials$year <- as.Date(as.character(brewing_materials$year),format="%Y")
brewing_materials$type <- as.factor(brewing_materials$type)

brewing_materials %>%
  group_by(year,type, month_current) %>%
  filter(str_detect(type, "Total" )) %>%
  summarize() %>%
  ggplot(aes(x=as.factor(type), y=month_current)) +
  geom_boxplot() 

eh <- brewing_materials %>%
  group_by(year,type, month_current) %>%
  filter(str_detect(type, "Total" )) %>%
  group_by(year, type) %>%
  summarise(summed = sum(month_current)) 

##this is weird why the fall after 2015?
b <- brewing_materials %>%
  group_by(year,type, month_current) %>%
  filter(str_detect(type, "Total" )) %>%
  filter(!year > "2015-03-31") %>%
  group_by(year, type) %>%
  summarise(summed = sum(month_current)) %>%
  ggplot() +
  geom_line(aes(x=year, y=summed, color=type), size=3)+
  scale_x_date(date_labels = "%Y")+
  theme_few() +
  labs(x = " ",
       y = "Raw Material In Pounds")+
  scale_y_continuous(label = ff_denom()) +
  scale_color_rickandmorty() + 
  theme(plot.title = element_text(size = 10))+
  ggtitle("I wonder why the grain products are being reduced ..?")


###

c <- brewing_materials %>%
  group_by(year, material_type, type) %>%
  filter(!str_detect(type, "Total")) %>%
  filter(!year > "2015-03-31") %>%
  summarise(dontknw = sum(month_current)) %>%
  filter(!str_detect(type, "Total")) %>%
  ggplot() +
  geom_line(aes(x=year, y=dontknw, linetype=material_type, color=type))+
  scale_x_date(date_labels = "%Y")+ 
  labs(x = " ",
       y = "Raw Material In Pounds")+
  theme_few()+ 
  scale_y_continuous(label = ff_denom()) +
  scale_color_rickandmorty() +
  ggtitle("Breakdown of the products") +
  theme(plot.title = element_text(size = 10))

#### 
# Produce image using graphics device

forBarPlot <- brewing_materials %>%
  group_by(year, material_type, type) %>%
  #filter(!year > "2015-03-31") %>%
  filter(!str_detect(type, 'Total')) %>%
  summarise(dontknw = sum(month_current)) %>%
  group_by(type, material_type) %>%
  filter(!str_detect(type, 'Total')) %>%
  summarise(total = sum(dontknw)) 

beerMug <- image_read("beer-mug.png")
frink2 <- image_scale(beerMug, "500")
frink2 <- frink2 %>%
  image_transparent(color="grey", fuzz=5)

#what does a beer contain on average 
d <- ggplot(forBarPlot, aes(1, total, fill=type)) +
  annotation_custom(rasterGrob(frink2, 
                               width = unit(0.8,"npc"),
                               height = unit(0.8,"npc")), 
                    -Inf, Inf, -Inf, Inf) +
  geom_bar(stat="identity", position="fill", alpha=0.9) +
  labs(x = " ", y = " Proportion of Ingredients") +
  scale_fill_rickandmorty()+
  scale_y_continuous(labels = ff_prop2percent(digits = 0))+
  theme_few()+
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        ) +
  ggtitle("What's in your average beer?") +
  theme(plot.title = element_text(size = 10))

dev.off()

#####


a <- na.omit(brewer_size) %>%
  mutate(used_up = total_shipped + taxable_removals,
         unused = total_barrels - used_up) %>%
  group_by(year) %>%
  summarise(totalUnused = sum(unused),
            totalUsed = sum(used_up)) %>%
  mutate(totalUnusedinGal = totalUnused*31,
         totalUsedinGal = totalUsed * 31) %>%
  select(year, totalUnusedinGal, totalUsedinGal) %>%
  melt(id.vars="year") %>%
  ggplot(aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  annotate("text", x = 1.5, 
           y = 7000000000,
           size=2,
           label = "Unused:
           Not accounted for by shipping or taxable numbers \n for that year")+
  labs(x= " ",
       y = "Beer in Gallons") + 
  scale_y_continuous(label = ff_denom()) +
  geom_curve(aes(x = 1.5, y = 5000000000, 
                 xend = 1, yend = 1000000000), 
             curvature = 0.4,
             arrow = arrow(length = unit(0.03, "npc"))) + 
  guides(fill=FALSE) + 
  ggtitle('Annually Sold / Shipped Beer 
          in the US ~ 11 B Gallons (2008-2014)')+
  theme_few() +
  theme(plot.title = element_text(size = 10))
  
dev.off()

library(patchwork)

(d|(b/c)) + plot_annotation(
    title = 'Some Summaries',
    subtitle = 'My First #TidyTuesday!',
    theme = theme(plot.title = element_text(size = 12))
  )
ggsave("ttwk14.png")
