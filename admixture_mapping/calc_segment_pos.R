args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
lai.dir <- args[2]
local.anc.dir <- args[3]

local.anc <- read.delim(paste0(local.anc.dir, "/local_ancestry_segments.txt"))[,c(1:3)]
chr.anc <- local.anc[local.anc$chr.nr == chr,]

chr.anc$end.snp.nr <- cumsum(chr.anc$nr.snps)
chr.anc$start.snp.nr <- c(1, chr.anc$end.snp.nr[-length(chr.anc$end.snp.nr)] + 1)

snps <- read.table(
  paste0(lai.dir, "/chr", chr, "/chr", chr, "_local_ancestry_snps.txt"))[,1]
chr.anc$begin.pos <- snps[chr.anc$start.snp.nr]
chr.anc$end.pos <- snps[chr.anc$end.snp.nr]
chr.anc$pos <- round(rowMeans(chr.anc[,c(dim(chr.anc)[2]-1,dim(chr.anc)[2])]),0)

write.table(chr.anc, 
            paste0(lai.dir, "/local_ancestry_annot_chr", chr, ".txt"),
            sep="\t", quote=F, row.names=F, col.names=T)