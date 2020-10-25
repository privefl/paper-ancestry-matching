df <- readRDS("data/matched-pop.rds")
table(df$continent)
#   Africa      Asia    Europe   North America    Oceania   South America
#    12244     11421    445233            2618       1663            1327

library(dplyr)
library(ggplot2)

# Make sure country names are matching
df2 <- df %>%
  mutate(country = forcats::fct_recode(
    country,
    "Democratic Republic of the Congo" = "Congo",
    "Republic of Congo" = "Congo",
    "Burkina Faso" = "Burkina",
    "Madeira Islands" = "Madeira",
    "United Arab Emirates" = "Emirates",
    "Serbia/Montenegro" = "Serbia",
    "China" = "Hong Kong"
  ))

world_map <- map_data("world") %>%
  mutate(region = forcats::fct_recode(
    region,
    "United Kingdom" = "UK",
    "South Georgia and the South Sandwich Islands" = "South Sandwich Islands",
    "South Georgia and the South Sandwich Islands" = "South Georgia",
    "Myanmar (Burma)" = "Myanmar",
    "Republic of Kosovo" = "Kosovo",
    "Antigua and Barbuda" = "Antigua",
    "Antigua and Barbuda" = "Barbuda",
    "Serbia/Montenegro" = "Serbia"))

all_countries <- unique(df2$country)
world_regions <- unique(world_map$region)
not_found <- all_countries[!all_countries %in% world_regions]
table(df2$country)[not_found]
sapply(not_found, grep, x = world_regions, value = TRUE) %>%
  setNames(not_found) %>%
  { .[lengths(.) != 0] }
sapply(world_regions, grep, x = not_found, value = TRUE) %>%
  setNames(world_regions) %>%
  { .[lengths(.) != 0] }


p <- ggplot({
  df2 %>%
    group_by(country) %>%
    summarise(unmatched = if (n() > 30) 100 * mean(is.na(infer)) else NA) %>%
    arrange(desc(unmatched)) %>%
    print(n = Inf) %>%
    right_join(world_map, by = c("country" = "region"))
}) +
  bigstatsr::theme_bigstatsr() +
  # coord_equal() +
  geom_polygon(aes(long, lat, group = group, fill = unmatched),
               color = "gray90", size = 0.1) +
  scale_fill_viridis_c(direction = -1, limits = c(0, 100)) +
  labs(x = "Longitude", y = "Latitude", fill = "% not matched") +
  theme(legend.position = "top")
p

# Locations of 1000G populations
# https://www.internationalgenome.org/data-portal/population
igsr_pop <- bigreadr::fread2("data/igsr_populations.tsv")

fam2 <- bigreadr::fread2("tmp-data/1000G_phase3_common_norel.fam2")
(all_pop <- unique(fam2[3:5]))
ind <- match(all_pop$Population, igsr_pop$`Population code`)
df <- data.frame(pop  = all_pop$Population,
                 desc = all_pop$`Population Description`,
                 lat  = igsr_pop$`Population latitude`[ind],
                 lon  = igsr_pop$`Population longitude`[ind])

Gujarat <- c(22.309425, 72.136230)
df[df$pop == "GIH", c("lat", "lon")] <- Gujarat
df[df$pop == "ITU", c("lat", "lon")] <- c(16.5, 79.5)
df[df$pop == "STU", c("lat", "lon")] <- c(6.927079,	79.861244)
# Amsterdam <- c(52.377956, 4.897070)
# df[df$pop == "CEU", c("lat", "lon")] <- Amsterdam
# Edinburgh <- c(55.953251, -3.188267)
# London <- c(51.509865, -0.118092)
# df[df$pop == "GBR", c("lat", "lon")] <- (Edinburgh + London) / 2
# Mexico <- c(19.432608, -99.133208)
# df[df$pop == "MXL", c("lat", "lon")] <- Mexico
# ...

p + geom_point(aes(lon, lat), data = df, color = "red", size = 1.5)
# ggsave("figures/map-unmatched.pdf", width = 11, height = 7)
