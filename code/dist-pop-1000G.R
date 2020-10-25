library(bigsnpr)
library(ggplot2)
theme_set(theme_bigstatsr())
NCORES <- 15

#### Compute PCA on 1000G ####

bedfile <- download_1000G("tmp-data")
(obj.bed <- bed(bedfile))
plink2 <- download_plink2("tmp-data")

ind.norel <- runonce::save_run({
  rel <- snp_plinkKINGQC(
    plink2.path = plink2,
    bedfile.in = bedfile,
    thr.king = 2^-4.5,
    make.bed = FALSE,
    ncores = NCORES,
    extra.options = "--memory 5000"
  )
  str(rel)
  which(!obj.bed$fam$sample.ID %in% rel$IID2)
}, file = "tmp-data/ind-norel-1000G.rds")

obj.svd <- runonce::save_run(
  bed_autoSVD(obj.bed, ind.row = ind.norel, k = 25, ncores = NCORES),
  file = "tmp-data/SVD-1000G.rds"
)
plot(obj.svd, type = "scores", scores = 14:25, coeff = 0.6)
PC <- predict(obj.svd)


#### Compute multiple distances between 1000G populations ####

fam2 <- bigreadr::fread2(sub_bed(bedfile, ".fam2"))[ind.norel, ]
pop_1000G <- fam2$Population
split_pop <- split(rows_along(fam2), pop_1000G)
npop <- length(split_pop)

## Fst

af_grp <- lapply(split_pop, function(ind) {
  bed_MAF(obj.bed, ind.norel[ind], attr(obj.svd, "subset"), ncores = NCORES)
})
dist_fst <- matrix(0, npop, npop)
colnames(dist_fst) <- rownames(dist_fst) <- names(split_pop)
for (i in 2:npop) {
  for (j in 1:(i - 1)) {
    fst <- snp_fst(af_grp[c(i, j)], overall = TRUE)
    dist_fst[i, j] <- dist_fst[j, i] <- print(fst)
  }
}
round(dist_fst, 3)
apply(round(100 * dist_fst, 2), 1, function(x) list(sort(x)[2:4]))

hc <- hclust(as.dist(dist_fst), method = "single")
pdf("figures/heatmap-Fst-1000G.pdf", width = 8, height = 8)
heatmap(dist_fst, symm = TRUE, Rowv = as.dendrogram(hc))
dev.off()

purrr::iwalk(list(AFR = 1:7, EAS = 8:12, SAS = 13:17, AMR = 18:21, EUR = 22:26), ~ {
  library(magrittr)
  dist_fst[hc$order, hc$order][, .x] %>%
    xtable::xtable(caption = "",
                   label = paste0("tab:fst-", .y),
                   align = paste(c("|l|", rep("c|", length(.x))), collapse = ""),
                   digits = 4) %>%
    print(caption.placement = "top",
          hline.after = c(-1, cumsum(c(0, 7, 5, 5, 4, 5))),
          file = paste0("tables/fst-values-", .y, ".tex"))
})


## Bhattacharyya distance

bhattacharyya <- function(m1, m2, s1, s2) {

  m_diff <- m2 - m1
  s_mean <- (s1 + s2) / 2

  res <- drop(crossprod(m_diff, solve(s_mean, m_diff))) / 8 +
    log(abs(det(s_mean)) / sqrt(abs(det(s1) * det(s2)))) / 2

  sqrt(res)
}

all_maha <- lapply(split_pop, function(ind) {
  bigutilsr::covrob_ogk(PC[ind, ])
})

dist_bhattacharyya <- matrix(0, npop, npop)
row.names(dist_bhattacharyya) <- colnames(dist_bhattacharyya) <- names(all_maha)
for (j in seq_len(npop)) {
  for (i in seq_len(j - 1)) {
    dist_bhattacharyya[i, j] <- dist_bhattacharyya[j, i] <- bhattacharyya(
      m1 = all_maha[[i]]$center, m2 = all_maha[[j]]$center,
      s1 = all_maha[[i]]$cov, s2 = all_maha[[j]]$cov)
  }
}
round(dist_bhattacharyya)

hc <- hclust(as.dist(dist_bhattacharyya), method = "single")
pdf("figures/heatmap-bhattacharyya-1000G.pdf", width = 8, height = 8)
heatmap(dist_bhattacharyya, symm = TRUE, Rowv = as.dendrogram(hc))
dev.off()

qplot(dist_bhattacharyya, dist_fst) +
  scale_color_viridis_c() +
  geom_smooth(color = "red", method = "lm", formula = "y ~ x + 0") +
  labs(x = "Bhattacharyya distance between clusters", y = "Fst")
# ggsave("figures/compare-Bhattacharyya-to-Fst.pdf", width = 8, height = 6)


## Distance between centers

all_centers <- bigutilsr::geometric_median(PC, by_grp = pop_1000G)
dist_centers <- dist(all_centers)

hc <- hclust(dist_centers, method = "single")
pdf("figures/heatmap-centers-1000G.pdf", width = 8, height = 8)
heatmap(as.matrix(dist_centers), symm = TRUE, Rowv = as.dendrogram(hc))
dev.off()

qplot(as.matrix(dist_centers)^2, dist_fst) +
  scale_color_viridis_c() +
  geom_smooth(color = "red", method = "lm", formula = "y ~ x + 0") +
  labs(x = "Squared distance between centers of clusters", y = "Fst")
# ggsave("figures/compare-Euclidean-to-Fst.pdf", width = 8, height = 6)

all_plots <- lapply(c(2, 3, 4, 8, 16, 25), function(K) {
  all_centers <- bigutilsr::geometric_median(PC[, 1:K], by_grp = pop_1000G)
  dist_centers <- dist(all_centers)
  qplot(as.matrix(dist_centers)^2, dist_fst) +
    scale_color_viridis_c() +
    theme_bigstatsr(0.7) +
    geom_smooth(color = "red", method = "lm", formula = "y ~ x + 0") +
    labs(x = "Squared distance between centers of clusters", y = "Fst") +
    ggtitle(paste("With", K, "PCs"))
})
plot_grid(plotlist = all_plots, scale = 0.95)
# ggsave("figures/compare-Euclidean-to-Fst2.pdf", width = 12, height = 6)


## Shortest distance

dist_closest <- matrix(0, npop, npop)
row.names(dist_closest) <- colnames(dist_closest) <- names(split_pop)
for (j in seq_len(npop)) {
  for (i in seq_len(j - 1)) {
    dist_closest[i, j] <- dist_closest[j, i] <-
      min(nabor::knn(PC[split_pop[[i]], ], PC[split_pop[[j]], ], k = 1)$nn.dists)
  }
}
round(dist_closest)

hc <- hclust(as.dist(dist_closest), method = "single")
pdf("figures/heatmap-closest-1000G.pdf", width = 8, height = 8)
heatmap(dist_closest, symm = TRUE, Rowv = as.dendrogram(hc))
dev.off()

qplot(dist_closest, dist_fst) +
  scale_color_viridis_c() +
  geom_smooth(color = "red", method = "lm", formula = "y ~ x + 0") +
  labs(x = "Shortest distance between clusters", y = "Fst")
# ggsave("figures/compare-closest-to-Fst.pdf", width = 8, height = 6)
