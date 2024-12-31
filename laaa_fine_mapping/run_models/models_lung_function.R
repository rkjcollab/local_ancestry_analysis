runLaaaModelSummary <- function(model.frame) {
  return (summary(lm(pheno ~ group + age + sex + ht + bmi + RFMIX_GW_AFR + allele_dose + afr_dose + allele_afr_dose, data=model.frame)))
}

runNullModel <- function(model.frame) {
  return (lm(pheno ~ group + age + sex + ht + bmi + RFMIX_GW_AFR, data=model.frame))
}

runLaaaModel <- function(model.frame) {
  return (lm(pheno ~ group + age + sex + ht + bmi + RFMIX_GW_AFR + allele_dose + afr_dose + allele_afr_dose, data=model.frame))
}

runAlleleModelSummary <- function(model.frame) {
  return (summary(lm(pheno ~ group + age + sex + ht + bmi + RFMIX_GW_AFR + allele_dose, data=model.frame)))
}
