# SDS 20250722

# Script used to prepare LAAA pheno files. Pulling all from the main phenotype
# file updated by Dr. Ortgea:
# Complete_SARP_CSGA_3-4-21_AA_cases_12yr_and_older_cleaned_5-9-2024.xlsx).

# Updated 20250702 to run LAAA stratified by sex on admixture mapping results
# stratified by sex.

# Phenotypes:

# “ALL MaxFVC”
# “ALL MaxFEV1FVC ratio”
# "ALL bFEV/FVC Ratio"
# "ALL MaxFEV1"

# Covariates:
# “ALL HT”
# “ALL AGE (csga-SARP12 default SARP3)”
# “ALL SEX”
# study group (SARP12+CSGA vs. SARP3)
# BMI
# RFMix_GW_AFR, calculated from RFMix run on unimputed data

# Setup ------------------------------------------------------------------------

library(readxl)
library(tidyverse)
library(conflicted)

setwd(paste0(Sys.getenv("RKJCOLLAB"), "/Collabs/ortega"))

# Using names as indicated in complete pheno file:
# Set desired phenotypes (one file will be made per phenotype)
pheno_list <- c("maxFVC" = "ALL MaxFVC",
                "maxFEV1_FVC" = "ALL MaxFEV1FVC ratio",
                "bFEV_FVC" = "ALL bFEV/FVC Ratio",
                "maxFEV1" = "ALL MaxFEV1")

# Set desired covariates (one file made for all covariates)
cov_list <- c("ht" = "ALL HT",
              "age" = "ALL AGE (csga-SARP12 default SARP3)",
              "sex" = "ALL SEX",
              "group" = "ALL COHORT Based on COMBINED Analysis",
              "bmi" = "BMI(csga-SARP12 default. Yellow SARP3)",
              "RFMIX_GW_AFR" = "RFMIX_GW_AFR")

# File IDs need to be edited below due to mix of different formats.

# Load data --------------------------------------------------------------------

# Load complete updated pheno (updated by SDS in update_complete_pheno.R)
pheno_full <- read_xlsx(paste0(
  "data/pheno/",
  "Complete_SARP_CSGA_3-4-21_AA_cases_12yr_and_older_cleaned_5-9-2024.xlsx"),
  na = "NA")

# List of IDs from RFMix (unimputed)
# Note: has a mix of formats - NWD & other
id_list <- read_delim(
  "data/rfmix_unimp/output_co/chr1/chr1_local_ancestry_samples.txt",
  col_names = c("fid", "iid"), delim = " ")

# Load global ancestry calculated by SDS from unimputed data (calculated by SDS)
anc_unimp <- read_delim("data/gwide_anc_unimp/merged_gwide.txt")

# List of IDs from RFMix (WGS) to filter files to
id_list_to_filt <- read_delim(
  "data/rfmix_wgs/output_o/chr8/chr8_local_ancestry_samples.txt",
  col_names = c("fid", "iid"), delim = " ")

# Make pheno -------------------------------------------------------------------

# Filter pheno to columns for covariates
# n_distinct(pheno_full$SARP_CSGA_genetics1_FID)   # 519, expect 519 / 521
pheno_full_filt <- pheno_full %>%
  dplyr::select(SARP_CSGA_genetics1_FID, TopMed_ID,
                all_of(pheno_list),
                any_of(cov_list)) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::mutate(SARP_CSGA_genetics1_FID = as.character(SARP_CSGA_genetics1_FID))

# Split by ID format
id_list_nwd <- id_list %>% dplyr::filter(str_detect(iid, "NWD"))
id_list_other <- id_list %>% dplyr::filter(!str_detect(iid, "NWD"))

# Make pheno file
pheno_nwd <- left_join(
  id_list_nwd, pheno_full_filt,
  by = c("iid" = "TopMed_ID")) %>%
  dplyr::mutate(TopMed_ID = iid)
pheno_other <- left_join(id_list_other, pheno_full_filt,
                         by = c("fid" = "SARP_CSGA_genetics1_FID"))

pheno <- full_join(pheno_nwd, pheno_other)

# Add global ancestry estimates from RFMix
pheno_anc_unimp <- inner_join(
  anc_unimp,
  pheno,
  by = c("IID" = "iid"))

# Update group variable to SARP12+CSGA vs. SARP3
table(pheno_anc_unimp$group)
# CSGA SARP1-2   SARP3 
# 116     289     107
pheno_anc_unimp_mod <- pheno_anc_unimp %>%
  dplyr::mutate(group_admix_map_lg = ifelse(group == "CSGA", "SARP1-2", group))
table(pheno_anc_unimp_mod$group_admix_map_lg)
# SARP1-2   SARP3 
# 405     107

# Recode group for PLINK
pheno_anc_unimp_mod_recode <- pheno_anc_unimp_mod %>%
  dplyr::mutate(group = case_when(
    group_admix_map_lg == "SARP1-2" ~ 0,
    group_admix_map_lg == "SARP3" ~ 1)) %>%
  dplyr::rename(FID = fid)

# Now, remove all NA values
pheno_anc_unimp_mod_recode_filt <- pheno_anc_unimp_mod_recode %>%
  dplyr::filter(if_all((names(cov_list)), ~ !is.na(.)))
  # dplyr::select(FID, IID, TopMed_ID, sex, age, ht, bmi, group, RFMIX_GW_AFR)
table(pheno_anc_unimp_mod_recode_filt$group) 
# 0   1 
# 396 107 

# Filter to just those IDs with WGS RFMix results
pheno_anc_unimp_mod_recode_filt_2 <- pheno_anc_unimp_mod_recode_filt %>%
  dplyr::filter(TopMed_ID %in% id_list_to_filt$fid)
table(pheno_anc_unimp_mod_recode_filt_2$group) 
#   0   1 
#   245 106 
n_distinct(pheno_anc_unimp_mod_recode_filt_2$IID)  # 351 / 351

# Make pheno files -------------------------------------------------------------

# For combined analysis
# for (pheno in names(pheno_list)) {
#   df <- pheno_anc_unimp_mod_recode_filt_2 %>%
#     dplyr::select(FID, IID, TopMed_ID, !!pheno, all_of(names(cov_list))) %>%
#     dplyr::filter(!is.na(!!sym(pheno)))
#   write_tsv(df, paste0(
#     "data/pheno/SARP123_CSGA_",
#     nrow(df), "_laaa_", pheno, "_pheno.txt"
#   ))
# }

# For sex-stratified analysis
cov_list <- grep("ALL SEX", cov_list, value = T, invert = T)
for (pheno in names(pheno_list)) {
  df <- pheno_anc_unimp_mod_recode_filt_2 %>%
    dplyr::filter(sex == 1) %>%
    dplyr::select(FID, IID, TopMed_ID, !!pheno, all_of(names(cov_list))) %>%
    dplyr::filter(!is.na(!!sym(pheno)))
  write_tsv(df, paste0(
    "data/pheno/strat_by_sex/SARP123_CSGA_",
    nrow(df), "_laaa_", pheno, "_pheno_male.txt"
  ))
}
for (pheno in names(pheno_list)) {
  df <- pheno_anc_unimp_mod_recode_filt_2 %>%
    dplyr::filter(sex == 2) %>%
    dplyr::select(FID, IID, TopMed_ID, !!pheno, all_of(names(cov_list))) %>%
    dplyr::filter(!is.na(!!sym(pheno)))
  write_tsv(df, paste0(
    "data/pheno/strat_by_sex/SARP123_CSGA_",
    nrow(df), "_laaa_", pheno, "_pheno_female.txt"
  ))
}
