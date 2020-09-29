args <- commandArgs(trailingOnly = TRUE)
r2.in.file <- args[1]
out.file <- args[2]
min.rsq <- as.numeric(args[3])

r2.frame <- read.delim(r2.in.file)
failed.r2.frame <- r2.frame[r2.frame$R2 < min.rsq,]
failed.snps <- paste(failed.r2.frame$CHROM, 
                     failed.r2.frame$POS,
                     failed.r2.frame$REF,
                     failed.r2.frame$ALT, sep=":")
write.table(failed.snps, out.file,  sep="\t", quote=F, row.names=F, col.names=F)


