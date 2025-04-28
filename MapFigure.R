library(rnaturalearth)
library(sf)
library(ggplot2)
install.packages("devtools")
devtools::install_github("ropensci/rnaturalearthhires")

#### BC/Washington Map

# Load Map data

canada <- ne_states(country = "Canada", returnclass = "sf")
usa <- ne_states(country = "United States of America", returnclass = "sf")

# Combine and crop
region_of_interest <- rbind(
  canada %>% filter(name == "British Columbia"),
  usa %>% filter(name == "Washington")
)

# crop more
cropped_region <- st_crop(region_of_interest, xmin = -128, xmax = -121, ymin = 47.5, ymax = 51)

points <- data.frame(
  name = c("Hood Head Farm", "Loughborough Inlet (Site 2)", "East West Bay", "Talbot", "Bamfield"),
  lat = c(47.884, 50.714, 50.305, 50.170, 48.834),
  long = c(-122.614, -125.446, -125.108, -124.879, -125.141)
)


# Plot 
map1 <- ggplot(data = cropped_region) +
  geom_sf(fill = "grey", color = "black") +
  geom_point(data = points, aes(x = long, y = lat), color = "red", size = 2) +
  geom_text(data = points, aes(x = long, y = lat, label = name), nudge_y = 0.08, size = 3) +
  annotate("point", x = -125.452, y =50.692, color = "red", size = 2)+
  annotate("text", x=-124.6, y = 50.69,
           label="Loughborough Inlet (Site 1)", size = 3)+
  annotate("point", x = -125.520, y =50.651, color = "red", size = 2)+
  annotate("text", x=-125.3, y = 50.58,
           label="Interfor", size = 3)+
  annotate("text", x= -123, y =50.3, label = "British Columbia", size = 5)+
  annotate("text", x= -121.9, y =48.8, label = "Washington", size = 5)+
  labs(x = NULL, y = NULL, caption = "Data: Natural Earth") +
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

map1

### Regional Map
# Load Canada 
canada2 <- ne_states(country = "Canada", returnclass = "sf")

# Filter for British Columbia
bc <- canada[canada$name == "British Columbia", ]

#Vancouver
vancouver_bbox <- st_bbox(c(xmin = -123.3, xmax = -123, ymin = 49.2, ymax = 49.43), crs = st_crs(bc))
vancouver <- st_crop(bc, vancouver_bbox)

points2 <- data.frame(
  name = c("Girl in Wetsuit", "Lighthouse Park", "Sandy Cove Park", "Third Beach"),
  lat = c(49.302, 49.335, 49.342, 49.304),
  long = c(-123.126, -123.263, -123.226, -123.157)
)


# Plot 
map2 <- ggplot(data = vancouver) +
  geom_sf(fill = "grey", color = "black") +
  geom_point(data = points2, aes(x = long, y = lat), color = "red", size = 3) +
  annotate("text", x = -123.263, y= 49.328, label = "Lighthouse Park", size = 3 )+
  annotate("text", x = -123.226, y= 49.350, label = "Sandy Cove Park", size = 3 )+
  annotate("text", x = -123.126, y= 49.295, label = "Girl In a Wetsuit", size = 3 )+
  annotate("text", x = -123.157, y= 49.315, label = "Third Beach", size = 3 )+
  annotate("text", x= -123.15, y = 49.4, label ="North Vancouver", size =5)+
  annotate("text", x= -123.15, y = 49.25, label ="Vancouver", size =5)+
  labs(x = NULL, y = NULL, caption = NULL) +
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank())
map2

## combine maps
library(cowplot)

map1 <- map1 + coord_sf()
map2 <- map2 + coord_sf()

mapfig <- plot_grid(map1, map2, ncol = 2, labels = c("A", "B"), rel_widths = c(1.8, 1))
mapfig

