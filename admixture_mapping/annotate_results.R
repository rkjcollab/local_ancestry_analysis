args <- commandArgs(trailingOnly = TRUE)

fam.file <- args[1]
rfmix.dir <- args[2]
model.pheno <- args[3]
assoc.file <- args[4]

# Limit to just those IDs run for each pheno's model, order to match PLINK
# input into RFMix so can extract African ancestry
model.ids <- read.delim(model.pheno)[,c(1,2)]
fam.ids <- read.table(fam.file, stringsAsFactors = F)[,c(1,2)]
names(fam.ids) <- c("FID", "IID")
fam.ids$ORDER <- 1:dim(fam.ids)[1]
merged.ids <- merge(fam.ids, model.ids)
merged.ids <- merged.ids[order(merged.ids$ORDER),]
fam.pos <- merged.ids$ORDER
hap.pos.1 <- fam.pos + (fam.pos-1)
hap.pos.2 <- fam.pos + (fam.pos)
hap.pos <- c(hap.pos.1, hap.pos.2)
hap.pos <- hap.pos[order(hap.pos)]

frame <- data.frame()
for (chr in 1:22) {
  results <- read.table(assoc.file, head=T, stringsAsFactors = F)
  results <- results[results$TEST == "ADD",]
  results <- results[results$CHR == chr, ]
  results$SEG.NR <- unlist(strsplit(results$SNP, ":"))[seq(2,dim(results)[1]*2,2)]
  results <- results[,c("CHR", "BETA", "P", "SEG.NR")]
  
  # Flip the sign because the "minor" allele would reflect European ancestry,
  # and we want to reflect the dose effect in terms of 0/1/2 copies of African
  # ancestry (PLINK association output always in terms of the minor allele).
  results$BETA <- results$BETA * (-1)
  
  # Read in RFMix calls, subtract 1 from all
  calls <- read.table(paste0(
    rfmix.dir, "/chr", chr, "/chr", chr, "_local_ancestry.0.Viterbi.txt")) - 1
  prop.afr.anc <- rowSums(calls[,hap.pos])/length(hap.pos)
  seg.nr <- 1:dim(calls)[1]
  seg.anc <- data.frame(prop.afr.anc, seg.nr)
  results <- merge(results, seg.anc, by.x="SEG.NR", by.y="seg.nr")  
  
  # Get additional info about segments from local_ancestry_annot files
  seg.info <- read.table(paste0(rfmix.dir, "/local_ancestry_annot_chr", chr, ".txt"), 
                         head=T, stringsAsFactors = F)
  results <- merge(results, seg.info, by.x="SEG.NR", by.y="seg.nr")
  
  results$begin.end.bp <- paste(results$begin.pos, results$end.pos, sep="-")
  results$nr.bp <- results$end.pos - results$begin.pos
  results$P <- format(results$P, scientific = T, digits = 3)
  results$prop.afr.anc <- round(results$prop.afr.anc, 2)
  results <- results[, c("CHR", "SEG.NR", "begin.end.bp",
                         "nr.bp", "nr.snps", "prop.afr.anc", "BETA", "P")]
  frame <- rbind(frame, results)
}

names(frame) <- toupper(names(frame))
frame <- frame[order(frame$CHR, as.numeric(frame$SEG.NR)),]
write.table(frame, paste0(assoc.file, "_annot.txt"), sep="\t", quote=F, row.names=F, col.names=T)
