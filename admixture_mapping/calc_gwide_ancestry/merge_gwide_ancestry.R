args <- commandArgs(trailingOnly = TRUE)
rfmix.results.dir <- args[1]
admix.out.dir <- args[2]
out.dir <- args[3]

rfmix <- read.delim(paste0(rfmix.results.dir, "/rfmix_gwide.txt"))
adm <- read.delim(paste0(admix.out.dir, "/admixture_gwide.txt"))
merged <- merge(rfmix, adm)
write.table(merged, paste0(out.dir, "/merged_gwide.txt"),
            sep="\t", quote=F, row.names=F, col.names=T)
