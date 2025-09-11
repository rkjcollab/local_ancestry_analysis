args <- commandArgs(trailingOnly = TRUE)

fam.file.name <- args[1]
rfmix.dir <- args[2]
out.file.name <- paste0(rfmix.dir, "/rfmix_gwide.txt")

# Get optional 3rd argument, name of SNPs per window file. Only need if RFMix
# not run in batches. Code automatically looks for batch 1 SNPs per window file.
if (length(args) == 3) {
  snps.per.win.file.name <- args[3] 
} else {
  snps.per.win.file.name <- NA
}

# Get values constant across all chromosomes
fam.file <- read.table(fam.file.name, header = F,
                       col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
nr.indiv <- nrow(fam.file)
nr.haplos <- nr.indiv*2
total.nr.snps <- 0
total.nr.afr.snps <- as.matrix(rep(0, nr.haplos))

for (chr in 1:22) {
  anc <- read.table(paste0(
    rfmix.dir, "/chr", chr, "/chr", chr, "_local_ancestry.0.Viterbi.txt"))
  # Get number of SNPs from SNP list
  nr.snps.check <- dim(read.delim(paste0(
    rfmix.dir, "/chr", chr, "/chr", chr, "_local_ancestry_snps.txt"),
    head=F, sep=" "))[1]
  
  # Check if SNPs per window file passed in (RFMix not run in batches)
  if (!is.na(snps.per.win.file.name)) {
    snps.per.win.file <- paste0(
      rfmix.dir, "/chr", chr, "/", snps.per.win.file.name)
  } else {
    # Set name for first batch
    snps.per.win.file <- paste0(
      rfmix.dir, "/chr", chr, "/chr", chr,
      "_local_ancestry_batch1.0.SNPsPerWindow.txt")
  }
   
  # Check if SNPs per window file exists
  if (file.exists(snps.per.win.file)) {
    print("SNPsPerWindow file found either found for batch 1 or passed in by")
    print("name. Using for calculation.")
    snps.per.win <- read.table(snps.per.win.file)
  } else {
    print("No SNPsPerWindow file found for batch 1 or passed in. If did not run")
    print("RFMix in batches, re-run with updated path for this file. If did run")
    print("RFMix in batches and with -o output, this is expected and file will")
    print("be made with 1 SNP per window for calculation.")
    snps.per.win <- data.frame(V1 = rep(1, nr.snps.check))
  }
  
  # Get number of SNPs from SNPs per window file
  nr.snps <- sum(snps.per.win$V1)
  
  # Make sure SNP numbers agree (for when SNPs per window file is loaded)
  if (nr.snps != nr.snps.check) {
    print(paste0(
      "ERROR! Number of SNPs does not check out for chromosome: ",
      chr,". Should be ", nr.snps.check, "; is ", nr.snps))
  }
  
  # Calculate number of AFR SNPs
  nr.afr.snps <- (as.matrix(t(anc))-1)%*%as.matrix(snps.per.win)
  total.nr.afr.snps <- total.nr.afr.snps + nr.afr.snps
  total.nr.snps <- total.nr.snps + nr.snps
}

# Get IID from fam file (col)
frame <- data.frame(IID = fam.file$IID)

# Complete frame output
frame$RFMIX_GW_AFR <- NA
k <- 0
for (i in seq(1,nr.haplos,2)) {
    k <- k + 1
    frame$RFMIX_GW_AFR[k] <-
      (total.nr.afr.snps[i] + total.nr.afr.snps[i+1])/(2*total.nr.snps)
}

write.table(frame, out.file.name, sep="\t", quote=F, row.names=F, col.names=T)
