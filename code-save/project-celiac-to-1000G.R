library(bigsnpr)
library(ggplot2)

# snp_plinkQC(
#   plink.path = download_plink("tmp-data"),
#   prefix.in = "~/Bureau/Dubois2010_data/FinnuncorrNLITUK3hap550",
#   geno = 0.01, mind = 0.01,
#   autosome.only = TRUE
# )

obj.bed <- bed("../Dubois2010_data/FinnuncorrNLITUK3hap550_QC.bed")
bed.ref <- bed(download_1000G("tmp-data"))
ind.norel <- readRDS("tmp-data/1008G-ind-norel.rds")

system.time(
  proj <- bed_projectPCA(bed.ref, obj.bed, ind.row.ref = ind.norel,
                         k = 25, ncores = nb_cores(),
                         join_by_pos = FALSE)
) # 4 min on laptop

PC.ref <- predict(proj$obj.svd.ref)
proj2 <- proj$OADP_proj

plot_grid(plotlist = lapply(1:9, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
    geom_point(aes(proj2[, k1], proj2[, k2]), color = "blue") +
    theme_bigstatsr(0.5) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2))
}), nrow = 3)

# ggsave("figures/proj-celiac-to-1000G.png", width = 11, height = 7)


#### Ancestry estimation ####

fam2 <- bigreadr::fread2("tmp-data/1000G_phase3_common_norel.fam2")[ind.norel, ]

source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
plot_grid2(plotlist = lapply(5:10, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2), color = fam2$Population) +
    theme_bigstatsr(0.7) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2))
}), title_ratio = 0, legend_ratio = 0.25)


seq_PC <- 1:18
pop <- paste(fam2$`Super Population`, fam2$Population, sep = "_")
pop_PCs <- vctrs::vec_split(PC.ref, pop)

all_pval <- sapply(pop_PCs$val, function(PC) {
  maha <- bigutilsr::covRob(PC[, seq_PC], estim = "pairwiseGK")
  dist <- mahalanobis(proj2[, seq_PC], center = maha$center, cov = maha$cov)
  pval <- pchisq(dist, df = length(seq_PC), lower.tail = FALSE)
})

choose_pop <- apply(all_pval, 1, which.max)
pval_max <- all_pval[cbind(rows_along(all_pval), choose_pop)]
choose_pop2 <- ifelse(pval_max < 0.05, NA, pop_PCs$key[choose_pop])

ggplot() +
  geom_histogram(aes(pval_max), breaks = seq(0, 1, by = 0.05),
                 color = "#FFFFFF", fill = "#000000", alpha = 0.5) +
  geom_vline(xintercept = 0.05, color = "red") +
  theme_bigstatsr() +
  labs(x = "Maximum p-value testing population belonging")

all_dist <- sapply(pop_PCs$val, function(PC) {
  maha <- bigutilsr::covrob_ogk(PC[, seq_PC])
  dist <- apply(proj2[, seq_PC], 1, function(x) sum((x - maha$center)^2))
})
choose_pop3 <- apply(all_dist, 1, which.min)
dist_min <- all_dist[cbind(seq_along(choose_pop3), choose_pop3)]
THR <- quantile(all_dist, 0.99) / 0.2 * 0.005
hist(log(dist_min)); abline(v = log(THR), col = "red")
choose_pop4 <- ifelse(dist_min > THR, NA, pop_PCs$key[choose_pop3])


# Get population from external files
pop.files <- list.files(path = "../Dubois2010_data/",
                        pattern = "cluster_*", full.names = TRUE)
pop <- snp_getSampleInfos(obj.bed, pop.files)[[1]]
pop_celiac <- c("Netherlands", "Italy", "UK", "UK", "Finland")[pop]

mean(is.na(choose_pop2)) # 29.6%
round(100 * tapply(is.na(choose_pop2), pop_celiac, mean), 1)
#   Finland       Italy Netherlands          UK
#      34.6        74.9        29.2        21.0

mean(is.na(choose_pop4)) # 0.09%

library(magrittr)
table(choose_pop4, pop_celiac, exclude = NULL) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Ancestry (left) of Celiac individuals and their matching to 1000G populations (top) by our method. See the description of 1000G populations at \\url{https://www.internationalgenome.org/category/population/}.",
                 label = "tab:ancestry-pred-celiac",
                 align = "|l|c|c|c|c|") %>%
  print(caption.placement = "top")

#         Finland Italy Netherlands   UK
# EUR_CEU      38     3        1161 4022
# EUR_FIN    2427     0           0    3
# EUR_GBR       2     0         472 2607
# EUR_IBS       1    89           5   34
# EUR_TSI       2   942          10   83
# <NA>          1     5           0    5

# \begin{table}[ht]
# \centering
# \caption{Ancestry (left) of Celiac individuals and their matching to 1000G populations (top) by our method. See the description of 1000G populations at \url{https://www.internationalgenome.org/category/population/}.}
# \label{tab:ancestry-pred-celiac}
# \begin{tabular}{|l|c|c|c|c|}
# \hline
# & Finland & Italy & Netherlands & UK \\
# \hline
# EUR\_CEU & 38 & 3 & 1161 & 4022 \\
# EUR\_FIN & 2427 &  &  & 3 \\
# EUR\_GBR & 2 &  & 472 & 2607 \\
# EUR\_IBS & 1 & 89 & 5 & 34 \\
# EUR\_TSI & 2 & 942 & 10 & 83 \\
# NA. & 1 & 5 &  & 5 \\
# \hline
# \end{tabular}
# \end{table}
