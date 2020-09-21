library(bigsnpr)
library(ggplot2)
bedfile <- download_1000G("tmp-data")
plink2 <- download_plink2("tmp-data")
rel <- snp_plinkKINGQC(
  plink2.path = plink2,
  bedfile.in = bedfile,
  thr.king = 2^-4.5,
  make.bed = FALSE,
  ncores = nb_cores()
)
str(rel)
(obj.bed <- bed(bedfile))
ind.rel <- match(c(rel$IID1, rel$IID2), obj.bed$fam$sample.ID)
ind.norel <- rows_along(obj.bed)[-ind.rel]
obj.svd <- bed_autoSVD(obj.bed, ind.row = ind.norel, k = 20,
                       ncores = nb_cores())
object.size(obj.bed)

fam2 <- bigreadr::fread2(sub_bed(bedfile, ".fam2"))

U <- obj.svd$u
pop_PCs <- vctrs::vec_split(U, with(fam2[ind.norel, ],
                                    paste0(`Super Population`, "_" ,Population)))
all_maha <- lapply(pop_PCs$val, bigutilsr::covrob_ogk)

bhattacharyya <- function(m1, m2, s1, s2) {

  m_diff <- m2 - m1
  s_mean <- (s1 + s2) / 2

  drop(crossprod(m_diff, solve(s_mean, m_diff))) / 8 +
    log(abs(det(s_mean)) / sqrt(abs(det(s1) * det(s2)))) / 2
}

npop <- length(all_maha)
dist <- matrix(0, npop, npop)
for (j in seq_len(npop)) {
  for (i in seq_len(j - 1)) {
    dist[i, j] <- dist[j, i] <- bhattacharyya(
      m1 = all_maha[[i]]$center, m2 = all_maha[[j]]$center,
      s1 = all_maha[[i]]$cov, s2 = all_maha[[j]]$cov
    )
  }
}
row.names(dist) <- colnames(dist) <- pop_PCs$key
dist

hc <- hclust(as.dist(sqrt(dist)))
plot(hc, xlab = "", ylab = "Distance", col = factor(substr(hc$labels, 1, 3), labels = viridis::viridis(5)))


dendrogram <- as.dendrogram(hc)
library(ggraph)
ggraph(dendrogram, 'dendrogram', circular = TRUE) +
  coord_equal() +
  geom_edge_elbow() +
  geom_node_text()

# install.packages("ape")
library("ape")
# Default plot
plot(as.phylo(hc), cex = 0.6, label.offset = 0.5)
library(dendextend)
labels_colors(dendrogram) <- 1
labels_colors(dendrogram) <- as.character(factor(substr(names(labels_colors(dendrogram)), 1, 3),
                                                 labels = scales::hue_pal()(5)))
labels_colors(dendrogram)
set(dendrogram, "branches_k_color", k = 2)
plot(
  set(dendrogram, "branches_k_color", k = 4), horiz = TRUE)

library("ggplot2")
library("ggdendro")
ggdendrogram(dendrogram, rotate = TRUE)


hcdata <- dendro_data(hc, type="rectangle")
hcdata$labels$superpop <- substr(hcdata$labels$label, 1, 3)

breaks <- c(sort(unique(round(segment(hcdata)$y, 1)), decreasing = TRUE)[1:4], 0)

ggplot() +
  theme_bw() +
  geom_segment(data=segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = label(hcdata), aes(x=x, y=y - 1, label=label, colour = superpop, hjust=0), size=3) +
  coord_flip() +
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position = "none") +
  scale_y_reverse(expand = c(0.2, 0), breaks = breaks, minor_breaks = seq(0, 160, by = 10)) +
  ggplot2::theme(
    rect = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor.y = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(y = "Distance")
