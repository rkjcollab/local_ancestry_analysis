args <- commandArgs(trailingOnly = TRUE)
in.file.name <- args[1]
chr <- args[2]

in.bim <- read.table(in.file.name, stringsAsFactors=F)
out.bim <- data.frame(V1=in.bim$V1,
                      V2=paste0(in.bim$V1, ":", in.bim$V4),
                      V3=rep(0,length(in.bim$V1)),
                      V4=in.bim$V4,
                      V5=in.bim$V5,
                      V6=in.bim$V6)

write.table(out.bim, paste0("tmp_new_chr", chr, ".bim"),  sep="\t", quote=F, row.names=F, col.names=F)
