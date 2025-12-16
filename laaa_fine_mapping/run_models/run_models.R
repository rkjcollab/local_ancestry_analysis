#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

pheno <- args[1]
chr <- args[2]
input.dir <- args[3]
output.dir <- args[4]
pheno.file <- args[5]
pheno_id_col_name <- args[6]
cov_list <- args[7]  # should be passed as "cov1,cov2"
code_dir <- args[8]

# Reformat covariate list
cov_list <- unlist(strsplit(cov_list, ","))

# Source lung function models
source(paste0(code_dir, "/run_models/models_lung_function.R"))

# Get the model frame & set indicated ID column to "id"
phenos <- read.delim(pheno.file, stringsAsFactors = F)
pheno.frame <- phenos
pheno.frame$id = pheno.frame[[pheno_id_col_name]]

# Load dose frames
allele.frame <- read.delim(
  paste0(input.dir, "/", pheno, "_chr_", chr, "_allele_dose.txt"), stringsAsFactors = F)
afr.frame <- read.delim(
  paste0(input.dir,  "/", pheno, "_chr_", chr, "_afr_dose.txt"), stringsAsFactors = F)
allele.afr.frame <- read.delim(
  paste0(input.dir,  "/", pheno, "_chr_", chr, "_allele_afr_dose.txt"), stringsAsFactors = F)

# Fix formatting to match pheno
colnames(allele.frame) <- gsub("\\.", ":", gsub("^X", "", colnames(allele.frame)))
colnames(afr.frame) <- gsub("\\.", ":", gsub("^X", "", colnames(afr.frame)))
colnames(allele.afr.frame) <- gsub("\\.", ":", gsub("^X", "", colnames(allele.afr.frame)))

# Create the output file
out.file.name <- paste0(output.dir, "/", pheno, "_assoc_chr", chr, ".txt")
cat(paste0(
  "position\tref\talt\talt_frq\tn\tunadj_allele_dose_beta\t",
  "unadj_allele_dose_beta_se\tunadj_allele_dose_p\tallele_dose_beta\t",
  "allele_dose_beta_se\tallele_dose_p\tafr_dose_beta\tafr_dose_beta_se\t",
  "afr_dose_p\tallele_afr_dose_beta\tallele_afr_dose_beta_se\t",
  "allele_afr_dose_p\tanova_p\n"),
  sep="", file=out.file.name, append=F)

# Loops once for each position in the allele frame
n <- NA
for (position in allele.frame$position) {
  # Merge in the allele dose
  # Makes new df converting dose frame ID columns to one column & allele dose
  # to the other column
  allele.col.frame <- 
    data.frame(id=names(allele.frame)[-c(1:3)],
              allele_dose=t(allele.frame[allele.frame$position == position,-c(1:3)]))
  names(allele.col.frame)[2] <- "allele_dose"
  model.frame <- merge(pheno.frame, allele.col.frame, by = "id")

  # Merge in the afr dose
  afr.col.frame <- 
    data.frame(id=names(afr.frame)[-c(1:3)],
              afr_dose=t(afr.frame[afr.frame$position == position,-c(1:3)]))
  names(afr.col.frame)[2] <- "afr_dose"
  model.frame <- merge(model.frame, afr.col.frame, by = "id")

  # Merge in the allele.afr dose
  allele.afr.col.frame <- 
    data.frame(id=names(allele.afr.frame)[-c(1:3)],
              allele.afr_dose=t(allele.afr.frame[allele.afr.frame$position == position,-c(1:3)]))
  names(allele.afr.col.frame)[2] <- "allele_afr_dose"
  model.frame <- merge(model.frame, allele.afr.col.frame, by = "id")

  # Update phenotype name
  names(model.frame)[names(model.frame) == pheno] <- "pheno"  

  # Remove variants not present (due to merging of SARP12 imputed and WGS)
  model.frame <- model.frame[!is.na(model.frame$allele_dose),]
  n <- dim(model.frame)[1]

  # Fit the model
  model <- runLaaaModelSummary(model.frame)
  m.null <- runNullModel(model.frame)
  m.laaa <- runLaaaModel(model.frame)
  anova_p <- anova(m.null, m.laaa)[2,6]
  allele.model <- runAlleleModelSummary(model.frame)

  # Get the ref and alt alleles
  ref <- allele.frame$ref[allele.frame$position == position]
  alt <- allele.frame$alt[allele.frame$position == position]

  # Estimate the alternate allele frequency
  frq <- sum(model.frame$allele_dose)/(dim(model.frame)[1]*2)

  # Get the model output
  if ("allele_dose" %in% rownames(model$coefficients)) {
    allele_dose_beta <- model$coefficients["allele_dose", "Estimate"]
    allele_dose_beta_se <- model$coefficients["allele_dose", "Std. Error"]
    allele_dose_p <- model$coefficients["allele_dose", "Pr(>|t|)"]      
  } else {
    allele_dose_beta <- NA
    allele_dose_beta_se <- NA
    allele_dose_p <- NA  
  }
  if ("afr_dose" %in% rownames(model$coefficients)) {
    afr_dose_beta <- model$coefficients["afr_dose", "Estimate"]
    afr_dose_beta_se <- model$coefficients["afr_dose", "Std. Error"]
    afr_dose_p <- model$coefficients["afr_dose", "Pr(>|t|)"]      
  } else {
    afr_dose_beta <- NA
    afr_dose_beta_se <- NA
    afr_dose_p <- NA  
  }
  if ("allele_afr_dose" %in% rownames(model$coefficients)) {
    allele_afr_dose_beta <- model$coefficients["allele_afr_dose", "Estimate"]
    allele_afr_dose_beta_se <- model$coefficients["allele_afr_dose", "Std. Error"]
    allele_afr_dose_p <- model$coefficients["allele_afr_dose", "Pr(>|t|)"]      
  } else {
    allele_afr_dose_beta <- NA
    allele_afr_dose_beta_se <- NA
    allele_afr_dose_p <- NA  
  }
  if ("allele_dose" %in% rownames(allele.model$coefficients)) {
    unadj_allele_dose_beta <- allele.model$coefficients["allele_dose", "Estimate"]
    unadj_allele_dose_beta_se <- allele.model$coefficients["allele_dose", "Std. Error"]
    unadj_allele_dose_p <- allele.model$coefficients["allele_dose", "Pr(>|t|)"]      
  } else {
    unadj_allele_dose_beta <- NA
    unadj_allele_dose_beta_se <- NA
    unadj_allele_dose_p <- NA  
  }

  # Write the output
  cat(position,"\t", ref, "\t", alt, "\t", frq, "\t", n, "\t",
      unadj_allele_dose_beta,"\t",
      unadj_allele_dose_beta_se,"\t",
      unadj_allele_dose_p,"\t",
      allele_dose_beta,"\t",
      allele_dose_beta_se,"\t",
      allele_dose_p,"\t",
      afr_dose_beta,"\t",
      afr_dose_beta_se,"\t",
      afr_dose_p,"\t",
      allele_afr_dose_beta,"\t",
      allele_afr_dose_beta_se,"\t",
      allele_afr_dose_p, "\t",
      anova_p, "\n", file=out.file.name, append=T, sep="")
}
