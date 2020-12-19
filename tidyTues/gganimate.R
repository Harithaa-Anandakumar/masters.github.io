library(dplyr)
library(tidyr)
library(tidyverse)
library(data.table)
library(reshape)
library(ggplot2)
library(forcats)
library(viridis)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(gghighlight)
library(plotly)
library(htmlwidgets)
library(gganimate)
library(stringr)
library(numform)

#read the data file 
annualDeaths <- read.csv("data/annual-number-of-deaths-by-cause-2.csv")

annualDeaths <- as.data.table(annualDeaths)
x <- data.table::melt(annualDeaths,
                      id.vars=c("Entity", "Code","Year"),
                      measure.vars= patterns("(deaths)"),
                      variable.name = "y"
)

fil <- x %>%
  dplyr::filter(Code %in% c("OWID_WRL"),
              #  Year %in% c("2017"),
                y != "Intestinal.infectious.diseases..deaths.")
fil <- fil %>%
  dplyr::filter(y != "Heat..hot.and.cold.exposure.deaths.")
fil$y <- gsub("..deaths.","\\", fil$y)

fil <- as.data.frame(fil)

fil %>%
  ggplot(aes(Year, value, group = y, color=y)) + 
  geom_line() + 
  geom_segment(aes(xend = 2017, yend = value, color=y), linetype = 2) + 
  geom_point(size = 2) + 
  geom_text(aes(x = 2017, label = y), hjust = 0) + 
  transition_reveal(Year) + 
  coord_cartesian(clip = 'off') + 
  scale_y_continuous(label = ff_denom()) +
  scale_colour_wsj() +
  labs(y = 'Number of Deaths in million', x = 'Year',
       title = 'Causes of death (1990-2016) in the world',
      # subtitle = 'Non-Communicable (~chronic), Communicable (~infectious+maternal+neonatal)',
       caption = "Data source: \n Global Burden of Disease Study 2017 (GBD 2017) Results.\n Downloaded data from OurWorldInData ") + 
  theme_minimal() + 
  guides(colour=FALSE)+
  theme(plot.margin = margin(1, 4, 2, 1, "cm"))

anim_save("00-animated-lineplot-Allcateg.gif")


