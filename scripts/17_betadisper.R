# ======================================================================
# 17_betadisper.R
#
# Homogeneity of multivariate dispersion among stress treatments: the check
# of the PERMANOVA equal-dispersion assumption described in the Methods.
#
# Uses exactly the inputs of the global PERMANOVA in
# scripts/02_fig1_diversity_permanova.R: the cleaned ASV relative-abundance
# matrix and aligned metadata from script 01 (eight-taxon community samples;
# ASVs with mean relative abundance <= 0.5% removed), the same Stress
# recoding, and the same Bray-Curtis distance matrix.
#
# Steps:
#   1. vegan::betadisper(): each sample's distance to its treatment centroid
#   2. vegan::permutest(), 999 permutations (as for the PERMANOVA in 02):
#      does mean distance to centroid differ among treatments?
#
# Output:
#   results/tables/Table_S29_betadisper.csv
#     one row per treatment: n, mean and SD of distance to centroid,
#     with the permutest F, df and P repeated on each row
#
# Run from the repository root:  Rscript scripts/17_betadisper.R
# Runtime: seconds. No model fitting.
# ======================================================================

source("scripts/utils_functions.R")
suppressPackageStartupMessages({ library(dplyr); library(vegan) })

ASV_t_clean        <- readRDS(P_RDS("ASV_matrix_clean.rds"))
env_option_A_clean <- readRDS(P_RDS("metadata_clean.rds"))

# Same Stress recoding as script 02
env_option_A_clean$Stress <- env_option_A_clean$Stress %>% as.character() %>%
  {replace(., . == "Salinity", "Sal")} %>% {replace(., . == "Temperature", "Temp")} %>%
  factor(levels = levels_short)

# Same distance matrix as script 02
Distances <- vegan::vegdist(ASV_t_clean, method = "bray")

# 1. Distance of each sample to its treatment centroid
mod <- vegan::betadisper(Distances, group = env_option_A_clean$Stress)

# 2. Permutation test of mean dispersion among treatments
set.seed(123)
pt <- vegan::permutest(mod, permutations = 999)
tab_pt <- pt$tab

F_stat <- tab_pt[["F"]][1]
df1    <- tab_pt[["Df"]][1]
df2    <- tab_pt[["Df"]][2]
p_val  <- tab_pt[["Pr(>F)"]][1]

res <- tibble::tibble(
  Stress   = as.character(env_option_A_clean$Stress),
  distance = as.numeric(mod$distances)
) %>%
  dplyr::group_by(Stress) %>%
  dplyr::summarise(n = dplyr::n(),
                   mean_distance_to_centroid = mean(distance),
                   sd_distance_to_centroid   = sd(distance),
                   .groups = "drop") %>%
  dplyr::mutate(Stress = factor(Stress, levels = levels_short)) %>%
  dplyr::arrange(Stress) %>%
  dplyr::mutate(Stress = as.character(Stress),
                mean_distance_to_centroid = round(mean_distance_to_centroid, 3),
                sd_distance_to_centroid   = round(sd_distance_to_centroid, 3),
                permutest_F   = round(F_stat, 3),
                permutest_df1 = df1,
                permutest_df2 = df2,
                permutest_P   = p_val,
                permutations  = 999)

readr::write_csv(res, P_TAB("Table_S29_betadisper.csv"))

cat(sprintf("permutest: F(%d,%d) = %.3f, P = %.3f (999 permutations)\n", df1, df2, F_stat, p_val))
print(as.data.frame(res[, c("Stress", "n", "mean_distance_to_centroid", "sd_distance_to_centroid")]))
message("Done: Table S29 written.")
