library(tidyverse)
library(ggrepel)
library(ggpubr)

# plot that lonely northern hemisphere

land <- map_data("world") %>% 
  filter(
    (region %in% c(
      "Greenland", 
      "Iceland", 
      "Canada", 
      "USA", 
      "Mexico", 
      "Guatemala",
      "Belize",
      "Honduras",
      "El Salvador",
      "Nicaragua",
      "Costa Rica"
      )
     ) |
      (subregion %in% c(
        "Svalbard"
        )
      )
  ) %>% 
  filter(long < 75)

ggplot() + 
  geom_polygon(
    data = land, 
    aes(x=long, y=lat, group=group),
    fill = "white"
  ) + 
  coord_map(projection = "ortho", orientation = c(45, -120, 0))
