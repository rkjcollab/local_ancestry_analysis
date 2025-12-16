args <- commandArgs(trailingOnly = TRUE)

bim.file.name  <- args[1] # Assumes that sample ID is in the 2nd column
rfmix.results.dir <- args[2]
admix.input.dir <- args[3]

bim <- read.table(bim.file.name, stringsAsFactors = F)

snps <- c()
for (chr in 22:22) {
# TEMP
# for (chr in 1:22) {
  positions <- read.table(paste0(
    rfmix.results.dir, "/chr", chr, "/chr", chr, "_local_ancestry_snps.txt"),
    head=F)[,1]
  chr.snps <- bim$V2[(bim$V1 == chr) & (bim$V4 %in% positions)]
  snps <- c(snps, chr.snps)
}

write.table(snps, paste0(admix.input.dir, "/initial_admix_snps.txt"),  
            sep="\t", quote=F, row.names=F, col.names=F)
