### Initialize document
# Load libraries
library(tidyverse)

# Load data
# soma_human <- readRDS("~/NNEDL/masterprojectxjdjg/curated/Human_studies/control-swap.rds") # includes both SomaHuman2 treatments
soma_human <- readRDS("~/NNEDL/masterprojectxjdjg/curated/Human_studies/control-swap-v2.rds") # only includes the large dose

################################################################################
#### Pre-process data
# Apply BH correction for multiple testing
soma_human$padj <- NA
for (norm in unique(soma_human$normalization_method)){
  for (analysis in unique(soma_human$analysis)){
    rows <- soma_human$normalization_method == norm & soma_human$analysis == analysis
    pval <- soma_human[rows,]$pvalue
    soma_human$padj[rows] <- p.adjust(pval, method = "BH")
  }
}

### Save data
save(soma_human,
     file = "~/NNEDL/masterprojectxjdjg/curated/Human_studies/soma_human.RData")


