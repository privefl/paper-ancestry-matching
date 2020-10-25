library(bigsnpr)
library(dplyr)
library(ggplot2)

print_table <- function(x) {
  tab <- sort(table(x, exclude = NULL), decreasing = TRUE)
  cat(paste0(names(tab), ": ", tab), sep = " || ")
}

#### Groupings of homogeneous ancestry  ####

df <- readRDS("data/info_UKBB.rds")
pop_UKBB <- df$pop
PC_UKBB <- as.matrix(select(df, PC1:PC16))

## 1) No known ancestry (just knowing there is an overdominant population)

round(one_center <- bigutilsr::geometric_median(na.omit(PC_UKBB)), 4)
#      PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8
# -11.6409   3.4545  -1.3283   0.7628  -1.2801  -0.2577   0.2160  -0.3044
#      PC9     PC10     PC11     PC12     PC13     PC14     PC15     PC16
#   0.2386   0.0030  -0.2657   0.1783   0.0167  -0.0744  -0.0025  -0.1004
round(bigutilsr::covrob_ogk(na.omit(PC_UKBB))$center, 4)
#      PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8
# -12.3377   3.7848  -1.5816   1.2128  -1.1087  -0.3453   0.2632  -0.3986
#      PC9     PC10     PC11     PC12     PC13     PC14     PC15     PC16
#   0.6842  -0.0162  -0.3409   0.2494   0.0271   0.0214  -0.0155  -0.1531

all_sq_dist <- bigutilsr:::rowSumsSq(sweep(PC_UKBB, 2, one_center, '-'))

qplot(log(all_sq_dist)) +
  theme_bigstatsr(0.9) +
  labs(x = "Log squared distance to overall center") +
  scale_x_continuous(breaks = 0:20)
# ggsave("figures/hist-dist-overall-center.pdf", width = 7, height = 5)

thr_sq_dist <- exp(7)  # visual inspection of histogram

mean(all_sq_dist < thr_sq_dist, na.rm = TRUE) # 91.06%
ind_homogeneous <- which(all_sq_dist < thr_sq_dist)

print_table(as.character(pop_UKBB[ind_homogeneous]))
# British: 421871 || Irish: 12039 || Other White: 8351 || NA: 1814 ||
# White: 467 || White and Asian: 23 || Indian: 5 ||
# White and Black African: 4 || White and Black Caribbean: 4 ||
# Black or Black British: 2 || African: 1 || Chinese: 1 || Other Asian: 1


## 2) Using self-reported ancestry

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

POP <- c("British", "Indian", "Pakistani", "Bangladeshi", "Chinese",
         "Caribbean", "African")
all_sq_dist <- apply(all_centers[POP, ], 1, function(center) {
  bigutilsr:::rowSumsSq(sweep(PC_UKBB, 2, center, '-'))
})
choose_pop <- apply(all_sq_dist, 1, function(x) {
  ind <- which.min(x)
  if (length(ind) == 0) NA else ind
})
min_sq_dist <- all_sq_dist[cbind(seq_along(choose_pop), choose_pop)]
(thr_sq_dist <- 0.002 * (max(dist(all_centers)^2) / 0.16))
hist(log10(min_sq_dist)); abline(v = log10(thr_sq_dist), col = "red")
choose_pop2 <- ifelse(min_sq_dist > thr_sq_dist, NA, POP[choose_pop])
mean(is.na(choose_pop2)) # 3.74%

table(pop_UKBB, factor(choose_pop2, levels = POP), exclude = NULL) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "",
                 label = "tab:ancestry-groups",
                 align = "|l|c|c|c|c|c|c|c|c|") %>%
  print(caption.placement = "top", hline.after = c(-1, 0, 4, 9, 12, 17, 18),
        sanitize.text.function = function(x) gsub("\\.", " ", x),
        file = "tables/ancestry-groups.tex")
#                           British Indian Pakistani Bangladeshi Chinese Caribbean African   <NA>
# British                    426210      6         4           1       1         2       0   4790
# Irish                       12712      0         0           0       0         0       0     41
# White                         492      0         0           0       1         0       0     52
# Other White                 10932      1         1           1       0         0       0   4880
# Indian                          6   1764      2488        1321       0         0       0    137
# Pakistani                       1    362      1299          63       0         0       0     23
# Bangladeshi                     0      3         0         215       0         0       0      3
# Chinese                         1      0         1           0    1437         0       0     65
# Other Asian                     4    113       169         745      62         0       1    653
# Caribbean                       0      2         0          23       0      2325    1148    799
# African                         1      0         1           0       0        74    2271    857
# Other Black                     0      1         1           1       0        36      33     46
# Asian or Asian British          0      7        16           3       1         0       0     15
# Black or Black British          2      0         0           0       0        11       9      4
# White and Black Caribbean       7      0         0           1       0        10       1    578
# White and Black African         6      0         0           0       0         1       2    393
# White and Asian                59     31         7          19       0         0       0    686
# <NA>                         2008    129       189         421     114       214     505   4240

df_pop0 <- df
df_pop0$pop <- choose_pop2
ind <- which(df_pop0$pop == "British")
df_pop <- rbind(df_pop0[-ind, ], df_pop0[sample(ind, 20e3), ])
alpha <- I(ifelse(is.na(df_pop$pop), 0.1, 1))

source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
theme_set(theme_bigstatsr(0.7))
plot_grid2(list(

  qplot(PC1, PC2, data = df_pop, color = pop, alpha = alpha) +
    labs(color = "Ancestry") +
    guides(colour = guide_legend(override.aes = list(size = 2))),

  qplot(PC3, PC4, data = df_pop, color = pop, alpha = alpha),
  qplot(PC5, PC6, data = df_pop, color = pop, alpha = alpha),
  qplot(PC7, PC8, data = df_pop, color = pop, alpha = alpha)

), scale = 0.95, nrow = 2, title_ratio = 0)
# ggsave("figures/UKBB-matched-ancestry.png", width = 9, height = 6)


## 3) Using country of birth

country2 <- forcats::fct_lump_min(df$country, 300)
print_table(country2)
# United Kingdom: 423572 || NA: 13865 || Ireland: 12525 || Other: 5850 ||
# India: 4012 || Caribbean: 2672 || Germany: 2136 || Kenya: 1684 ||
# Pakistan: 1439 || USA: 1390 || South Africa: 1364 || Nigeria: 1159 ||
# Ghana: 929 || Australia: 925 || France: 856 || Italy: 821 || Zimbabwe: 750 ||
# Sri Lanka: 744 || Canada: 729 || New Zealand: 686 || Malaysia: 679 ||
# Hong Kong: 648 || Poland: 637 || Uganda: 616 || Iran: 540 || Singapore: 502 ||
# Netherlands: 491 || Tanzania: 425 || China: 413 || Mauritius: 396 ||
# Malta: 365 || Spain: 355 || Iraq: 337 || The Guianas: 334 ||
# Philippines: 333 || Cyprus: 328 || Portugal: 320 || Egypt: 313 ||
# Barbados: 280 || Brazil: 271 || Japan: 266 || Zambia: 246 || Bangladesh: 246 ||
# Colombia: 245 || Denmark: 231 || Sierra Leone: 230 || Sweden: 216

all_centers <- bigutilsr::geometric_median(PC_UKBB, by_grp = country2)
dist_centers <- dist(all_centers)
attr(dist_centers, "Labels") <- glue::glue(
  "{rownames(all_centers)} ({table(country2)[rownames(all_centers)]})")
hc <- hclust(dist_centers, method = "single")
pdf("figures/heatmap-country-UKBB.pdf", width = 8, height = 8)
heatmap(as.matrix(dist_centers), symm = TRUE, Rowv = as.dendrogram(hc),
        margins = rep(8.5, 2))
dev.off()

# The clustering above guided the choice of these populations
POP <- c("United Kingdom", "Poland", "Iran", "Italy", "India",
         "China", "Caribbean", "Nigeria")
all_sq_dist <- apply(all_centers[POP, ], 1, function(center) {
  bigutilsr:::rowSumsSq(sweep(PC_UKBB, 2, center, '-'))
})
choose_pop <- apply(all_sq_dist, 1, function(x) {
  ind <- which.min(x)
  if (length(ind) == 0) NA else ind
})
min_sq_dist <- all_sq_dist[cbind(seq_along(choose_pop), choose_pop)]
(thr_sq_dist <- 0.002 * (max(dist(all_centers)^2) / 0.16))
hist(log10(min_sq_dist)); abline(v = log10(thr_sq_dist), col = "red")
choose_pop2 <- ifelse(min_sq_dist > thr_sq_dist, NA, POP[choose_pop])
mean(is.na(choose_pop2)) # 2.84%

table(pop_UKBB, factor(choose_pop2, levels = POP), exclude = NULL) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "",
                 label = "tab:country-groups",
                 align = "|l|c|c|c|c|c|c|c|c|c|") %>%
  print(caption.placement = "top", hline.after = c(-1, 0, 4, 9, 12, 17, 18),
        sanitize.text.function = function(x) gsub("\\.", " ", x),
        file = "tables/country-groups.tex")
#                           United Kingdom Poland   Iran  Italy  India  China Caribbean Nigeria   <NA>
# British                           423509   1412     30   3152     18      1         2       0   2890
# Irish                              12683     14      0     29      0      0         0       0     27
# White                                472     13      8     38      0      1         0       0     13
# Other White                         8102   2754    239   3259      2      0         0       0   1459
# Indian                                 6      0     33      0   4296      0         0       0   1381
# Pakistani                              1      0      2      0   1672      0         0       0     73
# Bangladeshi                            0      0      0      0      4      0         0       0    217
# Chinese                                1      0      0      0      0   1441         0       0     62
# Other Asian                            4      1    226      3    299     93         0       1   1120
# Caribbean                              0      0      0      0      3      0      2306    1245    743
# African                                1      0      0      0      2      0        71    2281    849
# Other Black                            0      0      0      0      2      0        36      34     46
# Asian or Asian British                 0      0      4      0     23      2         0       0     13
# Black or Black British                 2      0      0      0      0      0        11       9      4
# White and Black Caribbean              7      0      0      3      0      0        13       1    573
# White and Black African                6      0      0      4      0      0         1       2    389
# White and Asian                       56      0     12     30     54      0         0       0    650
# <NA>                                1827    116    680    462    345    315       215     513   3347

df_pop0 <- df
df_pop0$pop <- choose_pop2
ind <- which(df_pop0$pop == "United Kingdom")
df_pop <- rbind(df_pop0[-ind, ], df_pop0[sample(ind, 20e3), ])
alpha <- I(ifelse(is.na(df_pop$pop), 0.1, 1))

source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
theme_set(theme_bigstatsr(0.7))
plot_grid2(list(

  qplot(PC1, PC2, data = df_pop, color = pop, alpha = alpha) +
    labs(color = "Ancestry") +
    guides(colour = guide_legend(override.aes = list(size = 2))),

  qplot(PC3, PC4, data = df_pop, color = pop, alpha = alpha),
  qplot(PC5, PC6, data = df_pop, color = pop, alpha = alpha),
  qplot(PC7, PC8, data = df_pop, color = pop, alpha = alpha)

), scale = 0.95, nrow = 2, title_ratio = 0)
# ggsave("figures/UKBB-matched-country.png", width = 9, height = 6)

# Not matched
print_table(filter(df, PC4 < -80)$country)
# United Kingdom: 1460 || NA: 697 || USA: 150 || South Africa: 82 || Israel: 41 ||
# Egypt: 20 || Yemen: 17 || Hungary: 15 || Canada: 13
print_table(filter(df, PC6 > 50)$country)
# Colombia: 203 || Chile: 63 || Mexico: 57 || Peru: 53 || Ecuador: 35 || NA: 25 ||
# Venezuela: 23 || Bolivia: 22 || Brazil: 17 || Argentina: 10
print_table(filter(df, PC7 > 0 & PC8 < -30)$country)
# Colombia: 175 || Chile: 57 || Mexico: 51 || Peru: 51 || Ecuador: 33 ||
# Bolivia: 21 || NA: 16 || Venezuela: 15 || Brazil: 10
print_table(filter(df, PC7 < -20 & PC8 < -20)$country)
# United Kingdom: 363 || NA: 194 || USA: 28 || South Africa: 22 || Israel: 11
