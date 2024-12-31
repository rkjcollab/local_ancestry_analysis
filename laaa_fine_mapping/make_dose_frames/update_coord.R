args <- commandArgs(trailingOnly = TRUE)
bed <- args[1] # "$out_dir_hg38"/tmp_chr${chr}_out.bed
old.pheno.prefix <- args[2]  # "$out_dir_hg38"/${pheno}_chr_${chr}_
new.pheno.prefix <- args[3]
suffixes <- c("allele_dose.txt", "afr_dose.txt", "allele_afr_dose.txt")

bed.frame <- read.table(bed)[,c(2,4)]
names(bed.frame) <- c("new_position", "position")
for (suffix in suffixes) {
  file.old <- paste0(old.pheno.prefix, suffix)
  file.new <- paste0(new.pheno.prefix, suffix)
  dose.frame <- read.delim(file.old)
  dose.frame <- merge(dose.frame, bed.frame)
  dose.frame$position <- dose.frame$new_position
  dose.frame <- dose.frame[,-dim(dose.frame)[2]]
  write.table(dose.frame, file.new,  sep="\t", quote=F, row.names=F, col.names=T)
}