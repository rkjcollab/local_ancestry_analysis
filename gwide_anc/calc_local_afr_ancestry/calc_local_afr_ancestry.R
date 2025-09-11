args <- commandArgs(trailingOnly = TRUE)
nr.indiv <- as.integer(args[1])
nr.haplos <- nr.indiv*2
rfmix.results.dir <- args[2]
out.dir <- args[3]

frame <- data.frame()
for (chr in 1:22) {
   print(chr)
   calls <- read.table(paste0(
     rfmix.results.dir, "/chr", chr, "/chr", chr, "_local_ancestry.0.Viterbi.txt"))-1
   segments <- read.table(paste0(
     rfmix.results.dir, "/chr", chr, "/chr", chr, "_local_ancestry_batch1.0.SNPsPerWindow.txt"))
   prop.afr.anc <- rowSums(calls)/nr.haplos
   seg.nr <- 1:dim(segments)[1]
   chr.nr <- rep(chr, dim(segments)[1])
   nr.snps <- segments[,1]
   frame <- rbind(frame, data.frame(chr.nr, seg.nr, nr.snps, prop.afr.anc))
}

write.table(frame, paste0(
  out.dir, "/local_ancestry_segments.txt"),
  sep="\t", quote=F, row.names=F, col.names=T)
