library(bigsnpr)
library(ggplot2)
NCORES <- 15

#### Groupings of homogeneous ancestry  ####

## 1) Matching to 1000G populations

# Project UKBB individuals onto PCA space from 1000G

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

fam2 <- bigreadr::fread2("tmp-data/1000G_phase3_common_norel.fam2")[ind.norel, ]
pop_1000G <- paste0(fam2$`Super Population`, "_", fam2$Population)

seq_PC <- 1:16
thr_Fst <- 0.002

all_centers <- bigutilsr::geometric_median(PC.ref[, seq_PC], by_grp = pop_1000G)
thr_sq_dist <- thr_Fst * (max(dist(all_centers)^2) / 0.16)

all_sq_dist <- apply(all_centers, 1, function(center) {
  bigutilsr:::rowSumsSq(sweep(proj2[, seq_PC], 2, center, '-'))
})

## British ancestry (and around)
qplot(all_sq_dist[, "EUR_GBR"]) +
  theme_bigstatsr(0.9) +
  scale_x_log10() +
  labs(x = "Squared distance to GBR cluster") +
  geom_vline(xintercept = thr_sq_dist, color = "red")
# ggsave("figures/hist-dist-GBR.pdf", width = 7, height = 5)

mean(all_sq_dist[, "EUR_GBR"] < thr_sq_dist) # 89.73%
ind_GBR <- which(all_sq_dist[, "EUR_GBR"] < thr_sq_dist)

df0 <- readRDS("data/info_UKBB.rds")
df0 <- df0[match(obj.bed$fam$sample.ID, df0$eid), ]
pop_UKBB <- df0$pop
table(as.character(pop_UKBB[ind_GBR]), exclude = NULL)

## African ancestry
all_sq_dist_AFR <- all_sq_dist[, grepl("AFR", colnames(all_sq_dist))]
choose_pop_AFR <- apply(all_sq_dist_AFR, 1, which.min)
min_sq_dist_AFR <- all_sq_dist_AFR[cbind(seq_along(choose_pop_AFR), choose_pop_AFR)]
length(ind_AFR <- which(min_sq_dist_AFR < thr_sq_dist))  # 7443
table(as.character(pop_UKBB[ind_AFR]), exclude = NULL)

## East Asian
all_sq_dist_EAS <- all_sq_dist[, grepl("EAS", colnames(all_sq_dist))]
choose_pop_EAS <- apply(all_sq_dist_EAS, 1, which.min)
min_sq_dist_EAS <- all_sq_dist_EAS[cbind(seq_along(choose_pop_EAS), choose_pop_EAS)]
length(ind_EAS <- which(min_sq_dist_EAS < thr_sq_dist))  # 2314
table(as.character(pop_UKBB[ind_EAS]), exclude = NULL)

## South Asian
all_sq_dist_SAS <- all_sq_dist[, grepl("SAS", colnames(all_sq_dist))]
choose_pop_SAS <- apply(all_sq_dist_SAS, 1, which.min)
min_sq_dist_SAS <- all_sq_dist_SAS[cbind(seq_along(choose_pop_SAS), choose_pop_SAS)]
length(ind_SAS <- which(min_sq_dist_SAS < thr_sq_dist))  # 8213
table(as.character(pop_UKBB[ind_SAS]), exclude = NULL)

## Native American
all_sq_dist_AMR <- all_sq_dist[, grepl("AMR", colnames(all_sq_dist))]
choose_pop_AMR <- apply(all_sq_dist_AMR, 1, which.min)
min_sq_dist_AMR <- all_sq_dist_AMR[cbind(seq_along(choose_pop_AMR), choose_pop_AMR)]
length(ind_AMR <- which(min_sq_dist_AMR < thr_sq_dist))  # 216
table(as.character(pop_UKBB[ind_AMR]), exclude = NULL)

## Verification
ind_pop <- list(ind_GBR, ind_AFR, ind_EAS, ind_SAS, ind_AMR)
length(unlist(ind_pop)) # 456410
length(unique(unlist(ind_pop))) # 456410
matched_pop <- rep(c("GBR", "AFR", "EAS", "SAS", "AMR"), lengths(ind_pop))
df0_pop <- df0[unlist(ind_pop), ]

plot_grid(
  qplot(PC1, PC2, data = df0_pop, color = matched_pop) +
    theme_bigstatsr(0.8) +
    labs(color = "Ancestry") +
    theme(legend.position = c(0.75, 0.4)) +
    guides(colour = guide_legend(override.aes = list(size = 2))),

  qplot(PC3, PC4, data = df0_pop, color = matched_pop) +
    theme_bigstatsr(0.8) +
    theme(legend.position = "none"),

  scale = 0.95, nrow = 1
)

# ggsave("figures/UKBB-matched-ancestry.png", width = 11, height = 5)


## 2) Using self-reported ancestry (without projection)

library(dplyr)
PC_UKBB <- as.matrix(select(df0, PC1:PC16))
mixed <- c("Asian or Asian British", "Black or Black British",
           "White and Black Caribbean", "White and Black African",
           "White and Asian")
pop_UKBB2 <- forcats::fct_other(pop_UKBB, drop = mixed, other_level = NA) %>%
  droplevels(exclude = NA)
table(pop_UKBB2)
#     British       Irish       White Other White      Indian   Pakistani
#      431014       12753         545       15815        5716        1748
# Bangladeshi     Chinese Other Asian   Caribbean     African Other Black
#         221        1504        1747        4297        3204         118

all_centers <- bigutilsr::geometric_median(PC_UKBB, by_grp = pop_UKBB2)
thr_sq_dist <- thr_Fst * (max(dist(all_centers)^2) / 0.16)
dput(rownames(all_centers))
# c("British", "Irish", "White", "Other White", "Indian", "Pakistani",
#   "Bangladeshi", "Chinese", "Other Asian", "Caribbean", "African", "Other Black")

all_sq_dist <- apply(all_centers, 1, function(center) {
  bigutilsr:::rowSumsSq(sweep(PC_UKBB, 2, center, '-'))
})

## British ancestry (and around)
qplot(all_sq_dist[, "British"]) +
  theme_bigstatsr(0.9) +
  scale_x_log10() +
  labs(x = "Squared distance to British cluster") +
  geom_vline(xintercept = thr_sq_dist, color = "red")
# ggsave("figures/hist-dist-British.pdf", width = 7, height = 5)

mean(all_sq_dist[, "British"] < thr_sq_dist, na.rm = TRUE) # 92.66%
ind_British <- which(all_sq_dist[, "British"] < thr_sq_dist)
table(as.character(pop_UKBB[ind_British]), exclude = NULL)

## African ancestry
all_sq_dist_AFR <- all_sq_dist[, c("Caribbean", "African", "Other Black")]
length(ind_AFR2 <- which(
  apply(all_sq_dist_AFR, 1, function(x) any(x < thr_sq_dist))))  # 6767
table(as.character(pop_UKBB[ind_AFR2]), exclude = NULL)

## East Asian
length(ind_EAS2 <- which(all_sq_dist[, "Chinese"] < thr_sq_dist))  # 1616
table(as.character(pop_UKBB[ind_EAS2]), exclude = NULL)

## South Asian
all_sq_dist_SAS <- all_sq_dist[, c("Indian", "Pakistani", "Bangladeshi")]
length(ind_SAS2 <- which(
  apply(all_sq_dist_SAS, 1, function(x) any(x < thr_sq_dist))))  # 9409
table(as.character(pop_UKBB[ind_SAS2]), exclude = NULL)

## Native American
# No reported ancestry to match to


## 3) No known ancestry (just knowing there is an overdominant population)

one_center <- bigutilsr::geometric_median(na.omit(PC_UKBB))
all_sq_dist <- bigutilsr:::rowSumsSq(sweep(PC_UKBB, 2, one_center, '-'))

qplot(log(all_sq_dist)) +
  theme_bigstatsr(0.9) +
  labs(x = "Log squared distance to overall center") +
  scale_x_continuous(breaks = 0:20)
# ggsave("figures/hist-dist-overall-center.pdf", width = 7, height = 5)

thr_sq_dist <- exp(7)  # visual inspection of histogram

mean(all_sq_dist < thr_sq_dist, na.rm = TRUE) # 91.06%
ind_homogeneous <- which(all_sq_dist < thr_sq_dist)

table(as.character(pop_UKBB[ind_homogeneous]), exclude = NULL)


#### Compare indices ####

all_ind <- rlang::dots_list(ind_GBR, ind_British, ind_homogeneous,
                            ind_AFR, ind_AFR2, ind_EAS, ind_EAS2,
                            ind_SAS, ind_SAS2, ind_AMR, .named = TRUE)
names(all_ind) <- sub("ind_", "", names(all_ind))

lengths(all_ind)
#    GBR     British homogeneous         AFR        AFR2
# 438224      452441      444583        7443        6767
#    EAS        EAS2         SAS        SAS2         AMR
#   2314        1616        8213        9409         216

sapply(all_ind, function(ind) table(pop_UKBB[ind], exclude = NULL)) %>%
  rbind(t(data.frame(Total = lengths(all_ind)))) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "",
                 label = "tab:ancestry-groups",
                 align = "|l|c|c|c|c|c|c|c|c|c|c|") %>%
  print(caption.placement = "top", hline.after = c(-1, 0, 4, 9, 12, 17, 18, 19),
        sanitize.text.function = function(x) gsub("\\.", " ", x),
        file = "tables/ancestry-groups.tex")
#                              GBR British homogeneous  AFR AFR2  EAS EAS2  SAS SAS2 AMR
# British                   416374  426210      421871    2    2    1    1    6   11   0
# Irish                      12298   12712       12039    0    0    0    0    0    0   0
# White                        471     492         467    1    0    1    1    0    0   1
# Other White                 7187   10932        8351    0    0    0    0    1    3  40
# Indian                         5       6           5    0    0    0    0 4922 5573   0
# Pakistani                      0       1           0    0    0    0    0 1421 1724   0
# Bangladeshi                    0       0           0    0    0    0    0  217  218   0
# Chinese                        1       1           1    0    0 1453 1437    0    1   0
# Other Asian                    0       4           1    1    1  279   62  939 1027   0
# Caribbean                      0       0           0 3848 3578    0    0   25   25   0
# African                        1       1           1 2633 2347    0    0    0    1   0
# Other Black                    0       0           0   74   70    0    0    2    3   0
# Asian or Asian British         0       0           0    0    0    2    1   20   26   0
# Black or Black British         2       2           2   20   20    0    0    0    0   0
# White and Black Caribbean      7       7           4   24   14    0    0    1    1   1
# White and Black African        5       6           4    5    3    0    0    0    0   0
# White and Asian               21      59          23    0    0    2    0   26   57   1
# <NA>                        1852    2008        1814  835  732  576  114  633  739 173
# Total                     438224  452441      444583 7443 6767 2314 1616 8213 9409 216
