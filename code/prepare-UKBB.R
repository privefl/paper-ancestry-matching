##### Prepare genotypes ####

# file.symlink("~/NCRR-PRS/faststorage/UKBB/", ".")

write(sapply(1:22, function(chr) {
  paste0("UKBB/bed/",
         c(paste0("ukb_cal_chr", chr, "_v2.bed"),
           paste0("ukb_snp_chr", chr, "_v2.bim"),
           paste0("ukb58024_cal_chr", chr, "_v2_s488264.fam")))
}), tmp <- tempfile(), ncolumns = 3)

library(bigsnpr)
snp_plinkQC(
  plink.path = download_plink("tmp-data"),
  prefix.in = tmp,
  file.type = "--merge-list",
  prefix.out = "data/ukbb",
  geno = 0.01,
  autosome.only = TRUE,
  extra.options = "--memory 100000"
)

#### Prepare other data (self-reported ancestry and PCs) ####

library(bigreadr)
library(dplyr)

## Self-reported ancestry (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001)
code_ancestry <- fread2("UKBB/coding1001.tsv")
unknown <- c("Prefer not to answer", "Do not know", "Mixed",
             "Other ethnic group", "Any other mixed background")
code_country <- filter(fread2("UKBB/coding89.tsv"), selectable == "Y")
code_continent <- filter(fread2("UKBB/coding89.tsv"), selectable == "N")

df0 <- fread2(
  "UKBB/ukb41181.csv",
  select = c("eid", "21000-0.0", "20115-0.0", paste0("22009-0.", 1:16)),
  col.names = c("eid", "pop", "country", paste0("PC", 1:16))
) %>%
  mutate(
    pop  = factor(pop, levels = code_ancestry$coding,
                  labels = code_ancestry$meaning) %>%
      forcats::fct_recode("Other White" = "Any other white background",
                          "Other Asian" = "Any other Asian background",
                          "Other Black" = "Any other Black background") %>%
      forcats::fct_other(drop = unknown, other_level = NA) %>%
      droplevels(pop, exclude = NA) %>%
      forcats::fct_relevel(c(
        "British", "Irish", "White", "Other White",
        "Indian", "Pakistani", "Bangladeshi", "Chinese", "Other Asian",
        "Caribbean", "African", "Other Black")),

    continent = factor(country, levels = code_country$coding,
                       labels = code_country$parent_id) %>%
      factor(labels = code_continent$meaning),

    country = factor(country, levels = code_country$coding,
                     labels = code_country$meaning)
  )
str(df0)
df0$country[with(df0, is.na(country) & pop == "British")] <- "United Kingdom"
df0$country[with(df0, is.na(country) & pop == "Irish")]   <- "Ireland"
df0$continent[df0$country %in% c("United Kingdom", "Ireland")] <- "Europe"

bed_eid <- fread2("data/ukbb.fam")[[2]]
df <- df0[match(bed_eid, df0$eid), ]
mean(is.na(df$pop))  # 1.6%

saveRDS(df, "data/info_UKBB.rds")
