# SDS 20241211

# Set main paper colors.
library(colorspace)
library(khroma)

muted <- color("muted")

# Study colors
color_sarp <- muted(9)[2]  # "#332288"
color_csga <- muted(9)[3]  # "#332288"
color_both <- muted(9)[1]  # "#CC6677"

# Pheno colors
color_fvc <- base::unname(muted(9)[5])  # "#88CCEE"
color_fev1 <- base::unname(muted(9)[4])  # "#117733"
color_ratio <- base::unname(muted(9)[7])  # "#44AA99"

color_bfvc <- lighten(color_fvc, amount = 0.3)
color_bfev1 <- lighten(color_fev1, amount = 0.3)
color_bratio <- lighten(color_ratio, amount = 0.3)
color_maxfvc <- darken(color_fvc, amount = 0.3)
color_maxfev1 <- darken(color_fev1, amount = 0.3)
color_maxratio <- darken(color_ratio, amount = 0.3)
color_lfvc <- lighten(color_fvc, amount = 0.8)
color_lfev1 <- lighten(color_fev1, amount = 0.8)
color_lratio <- lighten(color_ratio, amount = 0.8)

# For SARP12 vs SARP3, want to use shades of SARP color
color_sarp12 <- lighten(color_sarp, amount = 0.3)  # "#6B62B2"
color_sarp3 <- darken(color_sarp, amount = 0.3)  # "#230085"