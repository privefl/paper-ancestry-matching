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

library(leaflet)
m <- leaflet(df) %>%
  addTiles() %>%
  addCircleMarkers(lat = ~lat, lng = ~lon, label = ~paste0(pop, ": ", desc),
                   radius = 2, color = "red") %>%
  print()

library(mapview)
mapshot(m, file = "figures/map-1000G.png")


coord <- cbind(df$lon, df$lat)
geo_dist <- sqrt(apply(coord, 1, geosphere::distGeo, p2 = coord))
rownames(geo_dist) <- colnames(geo_dist) <- df$pop
hc <- hclust(as.dist(geo_dist), method = "single")
heatmap(geo_dist, symm = TRUE, Rowv = as.dendrogram(hc))
