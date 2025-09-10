args <- commandArgs(trailingOnly = TRUE)

chr <- args[1]
work_dir <- args[2]
out_dir <- args[3]
perc_variance_exp <- as.numeric(args[4])

geno.fn <- paste0(work_dir, "/chr", chr, ".raw")
genos <- read.table(geno.fn, head=T)
genos <- genos[,7:dim(genos)[2]]
eigen.vals <- eigen(cor(genos))$values
denom <- sum(eigen.vals^2)
cumsums <- cumsum(eigen.vals^2/denom)
m.eff <- max(which(cumsums <= perc_variance_exp))
cat(chr, "\t", m.eff, "\t", dim(genos)[2], "\n",
    file=paste0(out_dir, "/m_eff_", perc_variance_exp, ".txt"), append = T)

