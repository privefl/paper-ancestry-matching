library(bigsnpr)
library(ggplot2)

snp_plinkQC(
  plink.path = download_plink("tmp-data"),
  prefix.in = "../POPRES_data/POPRES_allchr",
  geno = 0.01, mind = 0.01,
  autosome.only = TRUE
)
# map <- bigreadr::fread2("../POPRES_data/POPRES_allchr_QC.bim")
# tab <- bigreadr::fread2("../POPRES_data/POPRES_Snps_QC2.txt", header = FALSE)
# map$V2 <- tab$V2[match(map$V2, tab$V1)]
# bigreadr::fwrite2(map, "../POPRES_data/POPRES_allchr_QC.bim")

obj.bed <- bed("../POPRES_data/POPRES_allchr_QC.bed")
bed.ref <- bed(download_1000G("tmp-data"))
ind.norel <- readRDS("tmp-data/ind-norel-1000G.rds")

proj <- bed_projectPCA(bed.ref, obj.bed, ind.row.ref = ind.norel,
                       k = 25, ncores = nb_cores(),
                       join_by_pos = FALSE)
# 1,664,852 variants to be matched.
# 49,903 ambiguous SNPs have been removed.
# 275,628 variants have been matched; 49,165 were flipped and 87,081 were reversed.

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


#### Ancestry estimation with squared Euclidean distance on PCs ####

fam2 <- bigreadr::fread2("tmp-data/1000G_phase3_common_norel.fam2")[ind.norel, ]
pop_1000G <- paste0(fam2$`Super Population`, "_", fam2$Population)

plot_grid(plotlist = lapply(7:12, function(k) {
  k1 <- 2 * k
  k2 <- 2 * k + 1
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2), color = pop_1000G) +
    theme_bigstatsr(0.7) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    theme(legend.position = "none") +
    coord_equal()
}))
seq_PC <- 1:19

thr_Fst <- 0.002  # Fst between TSI and IBS is 0.0015

all_centers <- bigutilsr::geometric_median(PC.ref[, seq_PC], by_grp = pop_1000G)
thr_sq_dist <- thr_Fst * (max(dist(all_centers)^2) / 0.16)

all_sq_dist <- apply(all_centers, 1, function(center) {
  apply(proj2[, seq_PC], 1, function(x) sum((x - center)^2))
})
choose_pop <- apply(all_sq_dist, 1, which.min)
min_sq_dist <- all_sq_dist[cbind(seq_along(choose_pop), choose_pop)]
hist(log(min_sq_dist)); abline(v = log(thr_sq_dist), col = "red")
choose_pop2 <- ifelse(min_sq_dist > thr_sq_dist, NA, rownames(all_centers)[choose_pop])

pop <- obj.bed$fam$family.ID
pop2 <- dplyr::case_when(
  pop %in% c("Portugal", "Spain") ~ "SW Europe",
  pop %in% c("?orway", "Sweden", "Finland", "Denmark") ~ "Scandinavia",
  pop %in% c("Hungary", "Slovakia", "Czech", "Austria", "Slovenia", "Croatia") ~
    "Central Europe",
  pop %in% c("Russian", "Latvia", "Ukraine", "Poland") ~ "Eastern Europe",
  pop %in% c("Greece", "Turkey", "Serbia", "Cyprus", "Albania", "Kosovo",
             "Bosnia", "Macedonia,", "Romania", "Bulgaria") ~ "SE Europe",
  pop %in% c("Swiss-German", "Swiss-French", "Swiss-Italian") ~ "Switzerland",
  pop == "?etherlands" ~ "Netherlands",
  pop %in% c("United", "Scotland", "Ireland") ~ "Anglo-Irish Isles",
  TRUE ~ pop
)

mean(is.na(choose_pop2)) # 1.2%
round(100 * tapply(is.na(choose_pop2), pop2, mean), 1)
# Anglo-Irish Isles           Belgium    Central Europe    Eastern Europe
#               0.4               0.0               0.0               6.7
#            France           Germany             Italy       Netherlands
#               0.0               0.0               1.4               0.0
#       Scandinavia         SE Europe         SW Europe       Switzerland
#               0.0               9.6               0.4               0.0

library(magrittr)
table(pop2, choose_pop2, exclude = NULL) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Ancestry (left) of POPRES individuals and their matching to 1000G populations (top) by our method. See the description of 1000G populations at \\url{https://www.internationalgenome.org/category/population/}.",
                 label = "tab:ancestry-pred-popres",
                 align = "|l|c|c|c|c|c|c|") %>%
  print(caption.placement = "top", file = "tables/POPRES-matching.tex")

# https://www.internationalgenome.org/category/population/
#                   EUR_CEU EUR_FIN EUR_GBR EUR_IBS EUR_TSI <NA>
# Anglo-Irish Isles     136       0     127       0       2    1
# Belgium                43       0       0       0       0    0
# Central Europe         47       0       0       0       8    0
# Eastern Europe         27       0       0       0       1    2
# France                 49       0       3      35       2    0
# Germany                67       0       3       0       1    0
# Italy                   1       0       0      11     204    3
# Netherlands            13       0       4       0       0    0
# Scandinavia            13       1       1       0       0    0
# SE Europe              12       0       0       3      70    9
# SW Europe               1       0       0     261       1    1
# Switzerland           179       0       0      32      11    0
