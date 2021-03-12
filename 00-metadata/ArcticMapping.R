library(ggrepel)

latlon <- read_tsv("/Users/jwinnikoff/Downloads/arcticl_latlon.tsv")
arctic_coords <- latlon %>% 
  group_by(site) %>% 
  distinct(lat, lon, site)

# plot 'em

land <- map_data("world") %>% 
  filter(
    (region %in% c("Greenland", "Iceland")) |
      (subregion %in% c("Svalbard", "Alaska"))
  ) %>% 
  filter(long < 75)

arctic_coords %>% ggplot() + 
  geom_polygon(
    data = land, 
    aes(x=long, y=lat, group=group),
    fill = "white", colour = "black"
  ) + 
  geom_point(data = arctic_coords, aes(x=lon, y=lat), color="red") + 
  geom_text_repel(data = arctic_coords, aes(x=lon, y=lat, label=site)) +
  coord_map(projection = "globular")

# write for ODV
arctic_coords %>% 
  ungroup() %>% 
  dplyr::rename(LAT = lat, LONG = lon, `STATION` = site) %>% 
  mutate(
    CRUISE = "LINDBLAD2018",
    `Station ID` = row_number()
  ) %>% 
  write_csv("~/Downloads/ArcticCoords.csv")
