args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
info.file <- paste0("../../csga_imputation/data/output/chr", chr, ".info.gz")
out.file <- paste0("tmp_chr", chr, "_in.bed")

info <- read.delim(gzfile(info.file), stringsAsFactors = F)
info$Rsq_num <- as.numeric(info$Rsq)
info$CHR <- unlist(strsplit(info$SNP, split=":"))[seq(1,dim(info)[1]*4,4)]
info$POS <- as.numeric(unlist(strsplit(info$SNP, split=":"))[seq(2,dim(info)[1]*4,4)]) 
valid.snp.rows <- (((info$Genotyped == "Typed_Only") | (!is.na(info$Rsq_num) & (info$Rsq_num >= 0.7))) &
                  (info$MAF >= 0.01) & 
                  !(duplicated(info$POS))) 
out.frame <- info[,c("CHR", "POS", "POS", "SNP")]
write.table(out.frame, out.file, sep="\t", quote=F, row.names=F, col.names=F)

