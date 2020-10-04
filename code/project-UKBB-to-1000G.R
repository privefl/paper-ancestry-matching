library(bigsnpr)
library(ggplot2)
NCORES <- 15

#### Project UKBB individuals onto PCA space from 1000G ####

obj.bed <- bed("data/ukbb.bed")
bed.ref <- bed(download_1000G("tmp-data"))
ind.norel <- readRDS("tmp-data/ind-norel-1000G.rds")

proj <- runonce::save_run(
  bed_projectPCA(bed.ref, obj.bed, ind.row.ref = ind.norel,
                 k = 25, ncores = NCORES),
  file = "data/proj-UKBB-to-1000G.rds"
) # < 15 min
# 1,664,852 variants to be matched.
# 22,794 ambiguous SNPs have been removed.
# 377,914 variants have been matched; 0 were flipped and 65,346 were reversed.

PC.ref <- predict(proj$obj.svd.ref)
proj2 <- proj$OADP_proj

ind <- sample(nrow(proj2), 20e3)
# ind <- rows_along(proj2)

plot_grid(plotlist = lapply(1:9, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
    geom_point(aes(proj2[ind, k1], proj2[ind, k2]), color = "red", alpha = 0.3) +
    theme_bigstatsr(0.5) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    coord_equal()
}), nrow = 3)

# ggsave("figures/proj-UKBB-to-1000G.png", width = 8, height = 8)

fam2 <- bigreadr::fread2("tmp-data/1000G_phase3_common_norel.fam2")[ind.norel, ]
pop_1000G <- paste0(fam2$`Super Population`, "_", fam2$Population)
table(pop_1000G[PC.ref[, 9] > 40])  # 82 PUR
table(pop_1000G[PC.ref[, 14] > 60])  # 16 ITU


#### Ancestry estimation ####

plot(proj$obj.svd.ref)
plot_grid(plotlist = lapply(7:12, function(k) {
  k1 <- 2 * k
  k2 <- 2 * k + 1
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2), color = pop_1000G) +
    theme_bigstatsr(0.7) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    theme(legend.position = "none") +
    coord_equal()
}))
seq_PC <- 1:16

thr_Fst <- 0.002  # Fst between TSI and IBS is 0.0015

all_centers <- bigutilsr::geometric_median(PC.ref[, seq_PC], by_grp = pop_1000G)
thr_sq_dist <- thr_Fst * (max(dist(all_centers)^2) / 0.16)

all_sq_dist <- apply(all_centers, 1, function(center) {
  bigutilsr:::rowSumsSq(sweep(proj2[, seq_PC], 2, center, '-'))
})
choose_pop <- apply(all_sq_dist, 1, which.min)
min_sq_dist <- all_sq_dist[cbind(seq_along(choose_pop), choose_pop)]
hist(log(min_sq_dist)); abline(v = log(thr_sq_dist), col = "red")
choose_pop2 <- ifelse(min_sq_dist > thr_sq_dist, NA,
                      colnames(all_sq_dist)[choose_pop])
mean(is.na(choose_pop2)) # 4.55%


df0 <- readRDS("data/info_UKBB.rds")
pop_UKBB <- df0$pop[match(obj.bed$fam$sample.ID, df0$eid)]
mixed <- c("Asian or Asian British", "Black or Black British",
           "White and Black Caribbean", "White and Black African",
           "White and Asian")
pop_UKBB2 <- forcats::fct_other(pop_UKBB, drop = mixed, other_level = NA) %>%
  droplevels(exclude = NA)

dput(round(100 * tapply(is.na(choose_pop2), pop_UKBB, mean), 2))
# c(British = 2.22, Irish = 3.33, White = 7.89, `Other White` = 28.07,
#   Indian = 13.8, Pakistani = 18.71, Bangladeshi = 1.81, Chinese = 3.32,
#   `Other Asian` = 30.22, Caribbean = 9.87, African = 17.79, `Other Black` = 35.59,
#   `Asian or Asian British` = 47.62, `Black or Black British` = 15.38,
#   `White and Black Caribbean` = 94.3, `White and Black African` = 97.26,
#   `White and Asian` = 93.02)
dput(round(tapply(log10(min_sq_dist), pop_UKBB, mean), 1))
# c(British = 2.6, Irish = 2.7, White = 2.7, `Other White` = 2.9,
#   Indian = 2.7, Pakistani = 2.9, Bangladeshi = 2.1, Chinese = 2,
#   `Other Asian` = 2.7, Caribbean = 2.3, African = 2.3, `Other Black` = 2.8,
#   `Asian or Asian British` = 3.1, `Black or Black British` = 2.3,
#   `White and Black Caribbean` = 3.8, `White and Black African` = 3.8,
#   `White and Asian` = 3.6)

# https://www.internationalgenome.org/category/population/
library(magrittr)
table(pop_UKBB, substr(choose_pop2, 1, 3), exclude = NULL) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Self-reported ancestry (left) of UKBB individuals and their matching to 1000G continental populations (top) by our method. See the description of 1000G populations at \\url{https://www.internationalgenome.org/category/population/}.",
                 label = "tab:infer-UKBB-superpop",
                 align = "|l|c|c|c|c|c|c|") %>%
  print(caption.placement = "top", hline.after = c(-1, 0, 4, 9, 12, 17, 18),
        sanitize.text.function = function(x) gsub("\\.", " ", x),
        file = "tables/infer-UKBB-superpop.tex")
#                              AFR    AMR    EAS    EUR    SAS   <NA>
# British                        2      0      1 421457      6   9548
# Irish                          0      0      0  12328      0    425
# White                          1      1      1    499      0     43
# Other White                    0     40      0  11334      1   4440
# Indian                         0      0      0      5   4922    789
# Pakistani                      0      0      0      0   1421    327
# Bangladeshi                    0      0      0      0    217      4
# Chinese                        0      0   1453      1      0     50
# Other Asian                    1      0    279      0    939    528
# Caribbean                   3848      0      0      0     25    424
# African                     2633      0      0      1      0    570
# Other Black                   74      0      0      0      2     42
# Asian or Asian British         0      0      2      0     20     20
# Black or Black British        20      0      0      2      0      4
# White and Black Caribbean     24      1      0      8      1    563
# White and Black African        5      0      0      6      0    391
# White and Asian                0      1      2     27     26    746
# <NA>                         835    173    576   2296    633   3307

table(choose_pop2, pop_UKBB2, exclude = NULL) %>%
{ ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Self-reported ancestry (top) of UKBB individuals and their matching to 1000G populations (left) by our method. See the description of 1000G populations at \\url{https://www.internationalgenome.org/category/population/}.",
                 label = "tab:ancestry-fine-pred",
                 align = "|l|c|c|c|c|c|c|c|c|c|c|c|c|c|") %>%
  print(caption.placement = "top", hline.after = c(-1, 0, 7, 11, 16, 21, 26, 27),
        file = "tables/infer-UKBB-pop.tex")

## Verification if wrong matching
ind_mistmatch <- c(
  which(pop_UKBB %in% c("British", "Irish", "White") &
          substr(choose_pop2, 1, 3) != "EUR"),  # 12
  which(pop_UKBB %in% c("Indian", "Pakistani", "Bangladeshi",
                        "Chinese", "Other Asian") &
          substr(choose_pop2, 1, 3) != "SAS" &
          substr(choose_pop2, 1, 3) != "EAS"),  # 7
  which(pop_UKBB %in% c("Caribbean", "African", "Other Black") &
          substr(choose_pop2, 1, 3) != "AFR")  # 28, mainly Caribbean as SAS
) # 47 values
PC_UKBB <- df0[match(obj.bed$fam$sample.ID, df0$eid), -(1:2)]

ind <- sample(which(!is.na(pop_UKBB2)), 50e3)

theme_set(theme_bigstatsr(0.7))
PC_1_2 <- qplot(PC_UKBB[ind, 1], PC_UKBB[ind, 2],
                color = pop_UKBB2[ind], alpha = I(0.7)) +
  labs(x = "PC1", y = "PC2", color = "Self-reported ancestry") +
  scale_colour_hue(drop = FALSE)
PC_3_4 <- qplot(PC_UKBB[ind, 3], PC_UKBB[ind, 4],
                color = pop_UKBB2[ind], alpha = I(0.7)) +
  labs(x = "PC3", y = "PC4") +
  scale_colour_hue(drop = FALSE)

source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
plot_grid2(plotlist = list(

  PC_1_2,
  qplot(PC_UKBB[ind_mistmatch, 1], PC_UKBB[ind_mistmatch, 2],
        color = pop_UKBB2[ind_mistmatch], size = I(2), alpha = I(0.7)) +
    labs(x = "PC1", y = "PC2") +
    scale_x_continuous(limits = range(PC_UKBB[ind, 1], na.rm = TRUE)) +
    scale_y_continuous(limits = range(PC_UKBB[ind, 2], na.rm = TRUE)) +
    scale_colour_hue(drop = FALSE),

  PC_3_4,
  qplot(PC_UKBB[ind_mistmatch, 3], PC_UKBB[ind_mistmatch, 4],
        color = pop_UKBB2[ind_mistmatch], size = I(2), alpha = I(0.7)) +
    labs(x = "PC3", y = "PC4") +
    scale_x_continuous(limits = range(PC_UKBB[ind, 3], na.rm = TRUE)) +
    scale_y_continuous(limits = range(PC_UKBB[ind, 4], na.rm = TRUE)) +
    scale_colour_hue(drop = FALSE)

), title_ratio = 0, legend_ratio = 0.25, align = "hv")
# ggsave("figures/UKBB-mismatch.png", width = 10.5, height = 6.5)


#### Alternative ancestry matching: weighted 20-Nearest Neighbors ####
## Proposed in https://doi.org/10.1093/bioinformatics/btaa152

super_pop <- as.factor(substr(pop_1000G, 1, 3))

knn_dist <- bigutilsr::knn_parallel(PC.ref[, seq_PC], proj2[, seq_PC],
                                    k = 20, ncores = NCORES)
all_prob <- do.call("rbind", lapply(rows_along(proj2), function(i) {
  if (i %% 1000 == 0) print(i)
  w <- 1 / knn_dist$nn.dists[i, ]; w <- w / sum(w)
  prob <- tapply(w, super_pop[knn_dist$nn.idx[i, ]], sum)
}))
ind_ok <- which(all_prob > 0.875, arr.ind = TRUE)
pred_pop <- rep(NA, nrow(proj2))
pred_pop[ind_ok[, 1]] <- levels(super_pop)[ind_ok[, 2]]
table(pred_pop, exclude = NULL)

library(magrittr)
table(pop_UKBB, pred_pop, exclude = NULL) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Self-reported ancestry (left) of UKBB individuals and their matching to 1000G continental populations (top) using 20-wNN. See the description of 1000G populations at \\url{https://www.internationalgenome.org/category/population/}.",
                 label = "tab:ancestry-pred",
                 align = "|l|c|c|c|c|c|c|") %>%
  print(caption.placement = "top", hline.after = c(-1, 0, 4, 9, 12, 17, 18),
        sanitize.text.function = function(x) gsub("\\.", " ", x),
        file = "tables/matched-pop-KNN.tex")
#                              AFR    AMR    EAS    EUR    SAS   <NA>
# British                        4     50      6 430696     95    163
# Irish                          0      0      0  12748      3      2
# White                          1      2      1    540      1      0
# Other White                    0    170      1  15533     18     93
# Indian                         0      0      0     21   5680     15
# Pakistani                      0      0      0      3   1742      3
# Bangladeshi                    0      0      0      0    220      1
# Chinese                        0      7   1483      3      3      8
# Other Asian                    1      1    359    216   1138     32
# Caribbean                   4117      1      0      0     36    143
# African                     3000      0      1      2      2    199
# Other Black                   90      1      0      1      5     21
# Asian or Asian British         0      0      2      4     34      2
# Black or Black British        23      0      0      2      0      1
# White and Black Caribbean     93     16      0     74     11    403
# White and Black African      102     13      0     52      4    231
# White and Asian                0     42     10    242    349    159
# <NA>                        1024    541    712   3774   1020    749
