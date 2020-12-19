volcano <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-05-12/volcano.csv')
eruptions <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-05-12/eruptions.csv')
events <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-05-12/events.csv')
tree_rings <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-05-12/tree_rings.csv')
sulfur <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-05-12/sulfur.csv')
library(tidyverse)
library(dplyr)
library(tidyr)


sliced_eruptions  <- eruptions %>%
  arrange(desc(vei)) %>%
  select(volcano_number, vei, eruption_number, start_year, end_year) %>% 
  slice(1:100)


joined <- inner_join(volcano, sliced_eruptions, by="volcano_number")

joined <- inner_join(joined, events, by="volcano_number")
joined$volcano_name.x <- as.factor(joined$volcano_name.x)

joined <- joined %>%
  select(vei,volcano_name.x, primary_volcano_type,
         last_eruption_year, start_year, latitude,longitude,
         country, region, tectonic_settings,
         population_within_5_km:population_within_100_km) %>% 
  distinct()





library(maps)

world <- map_data("world", wrap=c(0,360))

library(wesanderson)
library(xkcd)

joined$lon2 <- ifelse(joined$longitude < -25, joined$longitude + 360, joined$longitude) # where d is your df
mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))




# map of volcano locations, colour coded by elevation

tiff("tidy20.tiff",units="in", width=8,height=5, res=300, compression = 'lzw') 
ggplot() +
  geom_polygon(
    data = mapWorld,
    aes(x = long, y = lat, group = group),
    fill = "grey",
    alpha = 0.4
  ) +
  geom_point(
    data = joined,
    aes(x = lon2, y = latitude, size=vei), shape=17,
    colour="darksalmon",
    alpha=0.6,
  ) +
  labs(title = "Top 100 destructive Volcanoes of All Time",
       subtitle = "Based on Volcanic Eruption Index-VEI",
       caption = "Most of the 'destructive' volcanoes are around the ring of fire. The Ring of Fire (also known as the Rim of Fire or the Circum-Pacific belt) 
       is a major area in the basin of the Pacific Ocean where many earthquakes and volcanic eruptions occur. In a large 40,000 km horseshoe shape, it is 
       associated with a nearly continuous series of oceanic trenches, volcanic arcs, and volcanic belts and plate movements. It has 452 volcanoes") +
  scale_size_continuous(breaks = c(5, 6, 7), labels = c("5","6","7"), name = "VEI")+
  theme(
    axis.line = element_blank(), 
    axis.title = element_blank(),
   # axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(colour = NA, fill="black"),
    panel.grid = element_blank(), 
   plot.subtitle = element_text(colour="mistyrose1", size=11),
   plot.caption = element_text(colour="mistyrose4", size=8),
   plot.title = element_text(colour="mistyrose"),
    plot.background = element_rect(fill = "black"),
    legend.title = element_text(colour='mistyrose2'),
    legend.background = element_rect(fill="black"),
    legend.key = element_rect(fill="black"),
    legend.text = element_text(colour="lightcoral"),
   text = element_text(family = "xkcd", colour = "lightcoral")
  )+
  annotate('text', x = 190, y = 20, 
           label = 'Ring of Fire', 
           color="mistyrose2",
           family="xkcd",
           size=6) +
  theme_xkcd()
dev.off()
