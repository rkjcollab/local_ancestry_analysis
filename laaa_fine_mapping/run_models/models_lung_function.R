runLaaaModelSummary <- function(model.frame) {
  form = as.character(paste(
    "pheno ~", paste(cov_list, collapse = " + "), "+ allele_dose + afr_dose + allele_afr_dose"))
  return (summary(lm(form, data=model.frame)))
}

runNullModel <- function(model.frame) {
  form = as.character(paste(
    "pheno ~", paste(cov_list, collapse = " + ")))
  return (lm(form, data=model.frame))
}

runLaaaModel <- function(model.frame) {
  form = as.character(paste(
    "pheno ~", paste(cov_list, collapse = " + "), "+ allele_dose + afr_dose + allele_afr_dose"))
  return (lm(form, data=model.frame))
}

runAlleleModelSummary <- function(model.frame) {
  form = as.character(paste(
    "pheno ~", paste(cov_list, collapse = " + "), "+ allele_dose"))
  return (summary(lm(form, data=model.frame)))
}
