args <- commandArgs(trailingOnly = TRUE)
in.file.name <- args[1]
chr <- args[2]

bim <- read.table(in.file.name, stringsAsFactors=F)
dupl.snps <- bim$V2[duplicated(bim$V2)]

write.table(dupl.snps, paste0("tmp_dupl_chr", chr, ".txt"),  sep="\t", quote=F, row.names=F, col.names=F)
