args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
begin.hg19 <- as.numeric(args[2])
end.hg19 <- as.numeric(args[3])
data_dir <- args[4]
out_dir <- args[5]

hg19pos <- read.delim(paste0(
  data_dir, "/snp_info.txt"),
  stringsAsFactors = F, sep=" ", head=F)[,1]

hg19.begin.index <- which(hg19pos == begin.hg19)
# Find the position before the specified beginning that is closest
if (length(hg19.begin.index) == 0) {
  hg19pos.delta <- begin.hg19 - hg19pos
  hg19pos.delta[hg19pos.delta < 0 ] <- NA
  hg19.begin.index <- which.min(hg19pos.delta)
  if (length(hg19.begin.index) == 0) {
    stop("Start of region for fine mapping not in input data.")
  }
}
hg19.end.index <- which(hg19pos == end.hg19)
# Find the position after the specified beginning that is closest
if (length(hg19.end.index) == 0) {
  hg19pos.delta <- hg19pos - end.hg19
  hg19pos.delta[hg19pos.delta < 0 ] <- NA
  hg19.end.index <- which.min(hg19pos.delta)
  if (length(hg19.end.index) == 0) {
    stop("End of region for fine mapping not in input data.")
  }
}

cat(hg19pos[hg19.begin.index],
    file=paste0(out_dir, "/tmp_", chr, "_begin.txt"), append=F)
cat(hg19pos[hg19.end.index],
    file=paste0(out_dir, "/tmp_", chr, "_end.txt"), append=F)

