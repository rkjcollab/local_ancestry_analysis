args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
file <- args[2]

bed <- read.table(file, stringsAsFactors=F)
bed <- bed[bed$V1 == paste0("chr", chr),]
bed <- bed[!duplicated(bed$V2),]
bed <- bed[order(bed$V2),]

write.table(bed, file,  sep="\t", quote=F, row.names=F, col.names=F)
