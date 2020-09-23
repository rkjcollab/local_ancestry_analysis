args <- commandArgs(trailingOnly = TRUE)
info.file <- args[1]
out.file <- args[2]
min.rsq <- as.numeric(args[3])
min.maf <- as.numeric(args[4])

info <- read.delim(gzfile(info.file), stringsAsFactors = F)
info$Rsq_num <- as.numeric(info$Rsq)
info$CHR <- unlist(strsplit(info$SNP, split=":"))[seq(1,dim(info)[1]*4,4)]
info$POS <- as.numeric(unlist(strsplit(info$SNP, split=":"))[seq(2,dim(info)[1]*4,4)]) 
valid.snp.rows <- (((info$Genotyped == "Typed_Only") | (!is.na(info$Rsq_num) & (info$Rsq_num >= min.rsq))) &
                  (info$MAF >= min.maf) & 
                  !(duplicated(info$POS))) 
out.frame <- info[,c("CHR", "POS", "POS", "SNP")]
write.table(out.frame, out.file, sep="\t", quote=F, row.names=F, col.names=F)

