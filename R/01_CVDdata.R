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

#read in data from "our world in data"
annualDeaths <- read.csv("data/annual-number-of-deaths-by-cause-2.csv")
totalDeaths <- read.csv("data/total-number-of-deaths-by-cause.csv")
totalDeaths <- as.data.table(totalDeaths)


#filter out only world data 
modtotalDeath <- totalDeaths %>%
  filter(Code == "OWID_WRL") 

colnames(modtotalDeath) <- c("Entity","Code","Year","Non-Communicable",
                 "Communicable","Injuries")
#gganimate plot
modtotalDeath %>%
  melt(id.vars=c("Entity","Code","Year")) %>%
  ggplot(aes(Year, value, group = variable, color=variable)) + 
  geom_line() + 
  geom_segment(aes(xend = 2017, yend = value, color=variable), linetype = 2) + 
  geom_point(size = 2) + 
  geom_text(aes(x = 2017, label = variable), hjust = 0) + 
  transition_reveal(Year) + 
  coord_cartesian(clip = 'off') + 
  scale_y_continuous(label = ff_denom()) +
  scale_colour_wsj() +
  labs(y = 'Number of Deaths in million', x = 'Year',
       title = 'Causes of death (1990-2016) in the world across 3 broad categories',
       subtitle = 'Non-Communicable (~chronic), Communicable (~infectious+maternal+neonatal)',
       caption = "Data source: \n Global Burden of Disease Study 2017 (GBD 2017) Results.\n Downloaded data from OurWorldInData ") + 
  theme_minimal() + 
  guides(colour=FALSE)+
  theme(plot.margin = margin(1, 3, 2, 1, "cm"))
anim_save("00-animated-lineplot-3categ.gif")

###still plot of the same
tiff("diseaseBreakdown.tiff",units="in", width=6,height=5, res=300, compression = 'lzw') 
modtotalDeath %>%
  melt(id.vars=c("Entity","Code","Year")) %>%
  ggplot(aes(Year, value, group = variable, color=variable)) + 
  geom_line() + 
  # geom_segment(aes(xend = 2017, yend = value, color=variable), linetype = 2) + 
  geom_point(size = 2) + 
  #geom_text(aes(x = 2017, label = variable), hjust = 0) + 
  #transition_reveal(Year) + 
  coord_cartesian(clip = 'off') + 
  scale_y_continuous(label = ff_denom()) +
  scale_colour_wsj() +
  labs(y = 'Number of Deaths in million', x = 'Year') +
  theme_few() +
  theme(legend.title = element_blank())

dev.off()

###still plot of the same ACCOUNTING for population by converting it to proportions
tiff("diseaseBreakdown.tiff",units="in", width=6,height=5, res=300, compression = 'lzw') 
modtotalDeath %>%
  melt(id.vars=c("Entity","Code","Year")) %>%
  ggplot(aes(Year, value, group = variable, color=variable)) + 
  geom_line() + 
  # geom_segment(aes(xend = 2017, yend = value, color=variable), linetype = 2) + 
  geom_point(size = 2) + 
  #geom_text(aes(x = 2017, label = variable), hjust = 0) + 
  #transition_reveal(Year) + 
  coord_cartesian(clip = 'off') + 
  scale_y_continuous(label = ff_denom()) +
  scale_colour_wsj() +
  labs(y = 'Number of Deaths in million', x = 'Year') +
  theme_few() +
  theme(legend.title = element_blank())

dev.off()

# Figure 1.2
tiff("diseaseBreakdown_02.tiff",units="in", width=6,height=5, res=300, compression = 'lzw') 
modtotalDeath %>%
  melt(id.vars=c("Entity","Code","Year")) %>%
  mutate(variable= factor(variable, levels=c("Communicable", "Non-Communicable","Injuries"))) %>% 
  ggplot(aes(Year, value, group = variable, fill=variable)) +
  geom_area(position = "fill", colour = "black", size = .2, alpha = .4) +
  #scale_fill_brewer(palette = "Blues") +
  scale_fill_wsj()+
  #geom_text(data=dd, aes(x=6, y=cum, label=name))
  scale_y_continuous(labels = scales::percent)+
  labs(y = 'Percentage proportion of deaths attributable by cause', x = 'Year') +
  theme_few() +
  theme(legend.title = element_blank())

dev.off()
annualDeaths <- as.data.table(annualDeaths)
x <- data.table::melt(annualDeaths,
                      id.vars=c("Entity", "Code","Year"),
                      measure.vars= pattern("(deaths)"),
                      variable.name = "y"
)

fil <- x %>%
  dplyr::filter(Code %in% c("OWID_WRL"),
                #Year %in% c("2017"),
                y != "Intestinal.infectious.diseases..deaths.")
fil$y <- gsub(".diseases","\\", fil$y)
fil <- fil %>%
  dplyr::filter(y != "Heat..hot.and.cold.exposure.")
tail(fil)
fil <- fil %>% mutate(rounded = value/10^6) %>%
  mutate_if(is.numeric, round, digits=3)

fil$y <- as.factor(fil$y)



#summing the total male & female suicides per country for each year
sm3<-aggregate(value~y+Year,fil,sum)
#* 1 ensures we have non-integer ranks while sliding
sm4<-sm3 %>% group_by(Year) %>% mutate(rank = min_rank(-value) * 1) %>%
  ungroup()

sm4 <- sm4 %>%
  mutate(rounded = (value/10^6)) %>%
  mutate_if(is.numeric, round)
sm4<-sm4 %>% group_by(Year) %>% mutate(rank = min_rank(-value) * 1) %>%
  ungroup()


static_plot <- ggplot(sm4,aes(rank,group=y,fill=as.factor(y),color=as.factor(y))) +
  geom_tile(aes(y = rounded/2,
                height = rounded,
                width = 0.9), alpha = 0.8, color = NA) +
  geom_text(aes(y = 0, label = paste(y, " ")), vjust = 0.2, hjust = 1) +
  geom_text(aes(y=rounded,label = paste(" ",rounded)), hjust=0)+
  coord_flip(clip = "off") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_reverse() +
  guides(color = FALSE, fill = FALSE) +
  theme_minimal() +
  theme(
    plot.title=element_text(size=12, hjust=0.5, face="bold", colour="grey", vjust=-1),
    plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="grey"),
    plot.caption =element_text(size=12, hjust=0.5, face="italic", color="grey"),
    axis.ticks.y = element_blank(), 
    axis.text.y = element_blank(), 
    plot.margin = margin(1, 1, 1, 5, "cm")
  )

library(gganimate)

plt <- static_plot + transition_states(
  states = Year, transition_length = 4, state_length = 1) + 
  ease_aes('cubic-in-out') +
 # view_follow(fixed_x = TRUE) +
  labs(title = "Total Deaths per Year per Disease: {closest_state}" 
       #subtitle = Top 10 Countries”,
       #caption = “Data Source: World Bank Data”,
       #x=”Countries”,y=”Total Suicides per year”)
  )

plt
animate(plt, nframes = 350,fps = 25,  width = 1200, height = 1000, 
        renderer = gifski_renderer("gganim.gif"))

#exploring more gganimate plots
#install.packages("gifs2ki")
library(gifski)
library(png)

animate(anim, nframes = 350,fps = 25,  width = 1200, height = 1000, 
        renderer = gifski_renderer("gganim.gif"))

anim_save("00-animated-barplot-transition.gif")

static_plot + transition_time(Year) +
  labs(title = "Look at the year: {closet_time}")


head(fil)
onlycar <- fil %>%
  filter(y == "Cardiovascular.diseases")
fill <- fil 
fil$value <- as.character(fil$value)
head(sm4)



ggplot(sm4, aes(x = y,y = value, fill = y)) +
  geom_bar(stat = "identity")  +
  geom_text(aes(label = value, 
                y = y),
            position = position_dodge(0.9),
            vjust = -1 ) +
  theme_classic() +
  transition_states(states = Year,
                    transition_length = 2, 
                    state_length = 1) + 
  labs(title = "Year: {closest_state}") +
  enter_fade() + 
  exit_shrink() +
  ease_aes('sine-in-out')
