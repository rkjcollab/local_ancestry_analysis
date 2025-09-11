args <- commandArgs(trailingOnly = TRUE)

fam.file.name <- args[1]
admix.fam.file.name <- args[2]
admix.in.dir <- args[3]
admix.out.dir <- args[4]
tgp.ids.file <- args[5]

# To note: ADMIXTURE output in same order as input file, have to use indexes

# Get PLINK study-only PLINK file and PLINK file that was input into ADMIXTURE
fam.file <- read.table(
  fam.file.name, header = F,
  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
admix.fam.file <- read.table(
  admix.fam.file.name, header = F,
  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))

# Get study IDs and number indv
nr.indiv <- nrow(fam.file)
aa.ids <- fam.file$IID
aa.rows <- which(!is.na(match(admix.fam.file$IID, aa.ids)))

# Get YRI rows matching pattern of IDs with current TGP input
tgp.ids <- read.delim(tgp.ids.file)
tgp.ids.afr <- tgp.ids$Sample.name[tgp.ids$Superpopulation.code == "AFR"]
yri.rows <- which(!is.na(match(admix.fam.file$IID, tgp.ids.afr)))
nr.yri <- sum(!is.na(yri.rows))

# Load ADMIXTURE results and select AFR column 
adm <- read.table(paste0(admix.out.dir, "/merged.2.Q"), head=F)
yri.anc <- adm[yri.rows, ]
if (mean(yri.anc$V1) > mean(yri.anc$V2)) {
  adm <- adm[, 1]
} else {
  adm <- adm[, 2]
}

# Add IDs from PLINK file
adm_ids <- data.frame(IID = admix.fam.file$IID, ADMIXTURE_GW_AFR = adm)

# Now, get ancestry just for study IDs
frame <- adm_ids[aa.rows, ]

write.table(frame, paste0(admix.out.dir, "/admixture_gwide.txt"),
            sep="\t", quote=F, row.names=F, col.names=T)
