#!/bin/bash

cd /Users/meherpreethi/Downloads/Victor_Ortega/SARP_CSGA_merge/pheno_analysis

/Users/meherpreethi/Downloads/plink_mac_20210606/plink --bfile ../rfmix/data/input/SARP123_CSGA_merged_final --allow-no-sex --covar SARP123_CSGA_512_covariates_num.txt --write-covar --dummy-coding --out SARP123_CSGA_512_covariates_recoded
