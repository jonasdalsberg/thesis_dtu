## Initialize document
library(tidyverse); library(magrittr); library(data.table); library(SomaDataIO)
source("scripts/SomaScan_Mouse/soma_norm_function.R")


# This script is to be loaded by the wrangle.R script in which the data is already loaded

################################################################################
## Define functions
#' Sample level median fold change normalization (alternative to ANML)
#'
#' @param ds data.table
#' @param v String of name of variable for values
#' @param ref String of name of variable for references
#' @param ext String giving extension given to scale factor variable name
#' @param mad String giving variable name of median absolute deviation of references
#' @param iteration Integer giving number of times to iterate
#' @importFrom magrittr %>%
#' @export
sample_norm <- function(ds, v, ref, ext) {
  sf_name <- paste0('sf', ext)
  ds[, (sf_name) := {
    cat('\rprogress', (.GRP / .NGRP) %>% multiply_by(100) %>% round(1))
    fifelse(SampleType %in% c('Sample', 'QC'),
            median(get(ref) / get(v)), 1)
  }, .(Dilution, SampleId, SlideId)]
}





################################################################################
## Cross-normalisation based on the recommendations from SomaLogic
unnorm1 <- read_adat("data/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.adat")
unnorm2 <- read_adat("data/raw/Mouse_studies/JYNR220601_GUS2022_898/SS-2231349_v4.1_other.adat")
hybNorm.medNormInt.plateScale_orig <- read_adat("data/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.hybNorm.medNormInt.plateScale.adat")
crossnorm <- soma_crossnorm(adat1 = unnorm1, adat2 = unnorm2, annotation_df = hybNorm.medNormInt.plateScale_orig)
crossnorm_JYNR210602 <- crossnorm[[1]]; crossnorm_JYNR220601 <- crossnorm[[2]]
# Align first df with wrangled data
crossnorm_JYNR210602 <- crossnorm_JYNR210602[crossnorm_JYNR210602$SampleType == "Sample",]
rownames(crossnorm_JYNR210602) <- crossnorm_JYNR210602$SampleId
crossnorm_JYNR210602 <- crossnorm_JYNR210602[rownames(crossnorm_JYNR210602) %in% rownames(log_JYNR210602),
                                             colnames(crossnorm_JYNR210602) %in% colnames(log_JYNR210602)]
crossnorm_JYNR210602 <- crossnorm_JYNR210602[order(rownames(crossnorm_JYNR210602)),]
sum(!rownames(crossnorm_JYNR210602) == rownames(log_JYNR210602)) # should be 0
crossnorm_JYNR210602 <- log2(crossnorm_JYNR210602)
# Align second df with wrangled data
crossnorm_JYNR220601 <- crossnorm_JYNR220601[crossnorm_JYNR220601$SampleType == "Sample",]
rownames(crossnorm_JYNR220601) <- crossnorm_JYNR220601$SampleId
crossnorm_JYNR220601 <- crossnorm_JYNR220601[rownames(crossnorm_JYNR220601) %in% rownames(log_JYNR220601),
                                             colnames(crossnorm_JYNR220601) %in% colnames(log_JYNR220601)]
crossnorm_JYNR220601 <- crossnorm_JYNR220601[order(rownames(crossnorm_JYNR220601)),]
sum(!rownames(crossnorm_JYNR220601) == rownames(log_JYNR220601)) # should be 0
crossnorm_JYNR220601 <- log2(crossnorm_JYNR220601)



################################################################################
## Cross-normalisation using MZCG's package (function)
# Load adats
ds1 <- read_adat("data/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.hybNorm.medNormInt.plateScale.adat")
ds2 <- read_adat("data/raw/Mouse_studies/JYNR220601_GUS2022_898/SS-2231349_v4.1_other.hybNorm.medNormInt.plateScale.adat")

# Extract SomaScan annotation metadata
meta <- ds1 %>% attributes() %>% .$Col.Meta %>% data.table()
meta[, variable := gsub("-", ".", SeqId) %>% paste0("seq.", .)]

# Get SOMAmer columns
seq_cols <- names(ds1) %>%  { .[grepl("^seq", .)] }

# Convert to data.tables
setDT(ds1); setDT(ds2)

# Put into long format
ds1.long <- ds1 %>% 
  melt(measure.vars = seq_cols) %>% 
  merge(meta, by = "variable")
ds2.long <- ds2 %>% 
  melt(measure.vars = seq_cols) %>%
  merge(meta, by = "variable")

# Derive reference values
ds1.long[, dataset := rep("JYNR210602"), .(variable)]
ds2.long[, dataset := rep("JYNR220601"), .(variable)]
ds.merged <- rbind(ds1.long, ds2.long)
ds.merged[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]

# Apply normalization
sample_norm(ds.merged, v = "value", ref = "ref.sample", ext = 5)
ds.merged[, value_norm := value * sf5]
ds.merged[, .(SampleId, variable, value, value_norm, dataset)] %>% 
  head()

# Split datasets again
ds1.long <- ds.merged[dataset == "JYNR210602" & SampleId %in% rownames(log_JYNR210602),]
ds2.long <- ds.merged[dataset == "JYNR220601" & SampleId %in% rownames(log_JYNR220601),]

# Sanity check (should be 0)
length(unique(ds1.long$SampleId)) - dim(log_JYNR210602)[1]
length(unique(ds2.long$SampleId)) - dim(log_JYNR220601)[1]

# Change to wide format and remove metadata
norm_JYNR210602 <- ds1.long[order(SampleId), .(SampleId, variable, value_norm)] %>% 
  dcast(SampleId ~ variable, value.var = "value_norm")
sum(!(norm_JYNR210602$SampleId == rownames(log_JYNR210602))) # should be 0
norm_JYNR210602 <- log2(as.data.frame(norm_JYNR210602)[,colnames(norm_JYNR210602) %in% colnames(log_JYNR210602)])
rownames(norm_JYNR210602) <- rownames(log_JYNR210602)

norm_JYNR220601 <- ds2.long[order(SampleId), .(SampleId, variable, value_norm)] %>% 
  dcast(SampleId ~ variable, value.var = "value_norm")
sum(!(norm_JYNR220601$SampleId == rownames(log_JYNR220601))) # should be 0
norm_JYNR220601 <- log2(as.data.frame(norm_JYNR220601)[,colnames(norm_JYNR220601) %in% colnames(log_JYNR220601)])
rownames(norm_JYNR220601) <- rownames(log_JYNR220601)





################################################################################
## Create within-study normalised data only based on cases or placebo
# (Makes for a more realistic scenario when cross-normalisation is to be applied,
#  as a real life example would not have a placebo group within study that is
#  almost identical to the one used in the cross-normalisation.
#  In a more realistic example, you would only have case (maybe several cases)
#  in the one study and then a placebo group (and possibly other cases) in the
#  second study.
#  Hence, we have the following potentially realistic scenarios:
  # Scenario 1:   Study 1: Case (DIO NASH vehicle), study 2: Placebo (Lean-CHOW vehicle)
  # Scenario 2:   Study 1: Cases (all DIO NASH groups), study 2: Placebo (Lean-CHOW vehicle)
  # Scenario 3:   Study 1: Case (DIO NASH vehicle), study 2: All groups
  # Scenario 4:   Study 1: Cases (all DIO NASH groups), study 2: All groups
# So we need to create and normalise three new dfs per study, i.e. 6 dfs in total
# And then we need to cross-normalise and then create 2 x 4 dfs (these are stored)
# For simplicity's sake, we can maybe make do with scenario 1 and scenario 4.
  # (Case/cases are going to be similar, as shown in the various volcano plots)
  # (So scenario 1 and 2 would be expected to produce somewhat similar results)
  # (Since DIO NASH vehicle is present in both studies, a cross-normalisation )
  # (like in scenario 3 seems like it would be slightly biased towards )
  # (producing too nice results.)

  # Also, new naming scheme.
  # Scenario 1 = scenario 1 (from above) with steps 1-4
  # Scenario 2 = scenario 4 (from above) with steps 1-4
  # Scenario 3 = scenario 1 (from above) with steps 1-5
  # Scenario 4 = scenario 4 (from above) with steps 1-5

# Procedure: Take 1-4 dfs and apply step 5 using MZCG's package...

# Put adats into long format
ds1.long <- ds1 %>% 
  melt(measure.vars = seq_cols) %>% 
  merge(meta, by = "variable")
ds2.long <- ds2 %>% 
  melt(measure.vars = seq_cols) %>%
  merge(meta, by = "variable")

# Save dataset information
ds1.long[, dataset := rep("JYNR210602"), .(variable)]
ds2.long[, dataset := rep("JYNR220601"), .(variable)]

# Filter for correct groups
ds1.long.case <- ds1.long[SampleType %in% c("Buffer", "Calibrator") | SampleId %in% rownames(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 3,]),]
ds1.long.all.cases <- ds1.long[SampleType %in% c("Buffer", "Calibrator") | SampleId %in% rownames(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group %in% c(3,4,5,6),]),]
ds1.long.placebo <- ds1.long[SampleType %in% c("Buffer", "Calibrator") | SampleId %in% rownames(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 1,]),]
ds1.long.nocont <- ds1.long[SampleType %in% c("Buffer", "Calibrator") | SampleId %in% rownames(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group %in% c(2,3,4,5,6),]),]
ds1.long.nonash <- ds1.long[SampleType %in% c("Buffer", "Calibrator") | SampleId %in% rownames(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group %in% c(1,2,4,5,6),]),]
ds2.long.case <- ds2.long[SampleType %in% c("Buffer", "Calibrator") | SampleId %in% rownames(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 2,]),]
ds2.long.all.cases <- ds2.long[SampleType %in% c("Buffer", "Calibrator") | SampleId %in% rownames(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group %in% c(2,3),]),]
ds2.long.placebo <- ds2.long[SampleType %in% c("Buffer", "Calibrator") | SampleId %in% rownames(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 1,]),]
ds2.long.nocont <- ds2.long[SampleType %in% c("Buffer", "Calibrator") | SampleId %in% rownames(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group %in% c(2,3),]),]
ds2.long.nonash <- ds2.long[SampleType %in% c("Buffer", "Calibrator") | SampleId %in% rownames(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group %in% c(1,3),]),]



# Derive reference values
ds1.long[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)] # all groups
ds1.long.case[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds1.long.all.cases[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds1.long.placebo[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds1.long.nocont[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds1.long.nonash[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds2.long[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds2.long.case[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds2.long.all.cases[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds2.long.placebo[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds2.long.nocont[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds2.long.nonash[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]


# Apply normalization
sample_norm(ds1.long, v = "value", ref = "ref.sample", ext = 5)
ds1.long[, value_norm := value * sf5]
sample_norm(ds1.long.case, v = "value", ref = "ref.sample", ext = 5)
ds1.long.case[, value_norm := value * sf5]
sample_norm(ds1.long.all.cases, v = "value", ref = "ref.sample", ext = 5)
ds1.long.all.cases[, value_norm := value * sf5]
sample_norm(ds1.long.placebo, v = "value", ref = "ref.sample", ext = 5)
ds1.long.placebo[, value_norm := value * sf5]
sample_norm(ds1.long.nocont, v = "value", ref = "ref.sample", ext = 5)
ds1.long.nocont[, value_norm := value * sf5]
sample_norm(ds1.long.nonash, v = "value", ref = "ref.sample", ext = 5)
ds1.long.nonash[, value_norm := value * sf5]
sample_norm(ds2.long, v = "value", ref = "ref.sample", ext = 5)
ds2.long[, value_norm := value * sf5]
sample_norm(ds2.long.case, v = "value", ref = "ref.sample", ext = 5)
ds2.long.case[, value_norm := value * sf5]
sample_norm(ds2.long.all.cases, v = "value", ref = "ref.sample", ext = 5)
ds2.long.all.cases[, value_norm := value * sf5]
sample_norm(ds2.long.placebo, v = "value", ref = "ref.sample", ext = 5)
ds2.long.placebo[, value_norm := value * sf5]
sample_norm(ds2.long.nocont, v = "value", ref = "ref.sample", ext = 5)
ds2.long.nocont[, value_norm := value * sf5]
sample_norm(ds2.long.nonash, v = "value", ref = "ref.sample", ext = 5)
ds2.long.nonash[, value_norm := value * sf5]




################################################################################
## Cross-normalisation using MZCG's package (function), based on only case and placebo
# Scenario 1:   Study 1: Case (DIO NASH vehicle), study 2: Placebo (Lean-CHOW vehicle), steps 1-4
# Scenario 2:   Study 1: Cases (all DIO NASH groups), study 2: All groups, steps 1-4
# Scenario 3:   Study 1: Case (DIO NASH vehicle), study 2: Placebo (Lean-CHOW vehicle), steps 1-5
# Scenario 4:   Study 1: Cases (all DIO NASH groups), study 2: All groups, steps 1-5
# Scenario 5:   Study 1: All groups except Lean-CHOW vehicle, study 2: All groups except DIO NASH vehicle, steps 1-4
# Scenario 6:   Study 1: All groups except Lean-CHOW vehicle, study 2: All groups except DIO NASH vehicle, steps 1-5

## JYNR210602, SCENARIO 1
  # Derive reference values and merge
  ds1.sce1 <- rbind(ds1.long.case, ds2.long.placebo)
  ds1.sce1[SampleType %in% c('Sample', 'QC'), ref.sample.cross := median(value), .(variable)]
  
  # Apply normalization
  sample_norm(ds1.sce1, v = "value", ref = "ref.sample.cross", ext = 6)
  ds1.sce1[, value_cross_norm := value * sf6]
  ds1.sce1[, .(SampleId, variable, value, value_norm, value_cross_norm, dataset)] %>% 
    head()
  
  # Split datasets again
  JYNR210602.sce1.norm.case <- ds1.sce1[dataset == "JYNR210602" & SampleId %in% rownames(log_JYNR210602),]
  JYNR220601.sce1.norm.placebo <- ds1.sce1[dataset == "JYNR220601" & SampleId %in% rownames(log_JYNR220601),]
  
  # Sanity check (should be 0)
  length(unique(JYNR210602.sce1.norm.case$SampleId)) - dim(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 3,])[1] # should be 0
  length(unique(JYNR220601.sce1.norm.placebo$SampleId)) - dim(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 1,])[1] # should be 0
  
  # Change to wide format and remove metadata
  JYNR210602.sce1.norm.case <- JYNR210602.sce1.norm.case[order(SampleId), .(SampleId, variable, value_cross_norm)] %>% 
    dcast(SampleId ~ variable, value.var = "value_cross_norm")
  JYNR210602.sce1.norm.case <- log2(as.data.frame(JYNR210602.sce1.norm.case)[,colnames(JYNR210602.sce1.norm.case) %in% colnames(log_JYNR210602)])
  rownames(JYNR210602.sce1.norm.case) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,])
  meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(JYNR210602.sce1.norm.case)]
  
  JYNR220601.sce1.norm.placebo <- JYNR220601.sce1.norm.placebo[order(SampleId), .(SampleId, variable, value_cross_norm)] %>% 
    dcast(SampleId ~ variable, value.var = "value_cross_norm")
  JYNR220601.sce1.norm.placebo <- log2(as.data.frame(JYNR220601.sce1.norm.placebo)[,colnames(JYNR220601.sce1.norm.placebo) %in% colnames(log_JYNR220601)])
  rownames(JYNR220601.sce1.norm.placebo) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,])
  meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(JYNR220601.sce1.norm.placebo)]


  
# ## JYNR210602, SCENARIO 2
#   # Derive reference values and merge
#   ds1.sce4 <- rbind(ds1.long.all.cases, ds2.long)
#   ds1.sce4[SampleType %in% c('Sample', 'QC'), ref.sample.cross := median(value), .(variable)]
# 
#   # Apply normalization
#   sample_norm(ds1.sce4, v = "value", ref = "ref.sample.cross", ext = 6)
#   ds1.sce4[, value_cross_norm := value * sf6]
#   ds1.sce4[, .(SampleId, variable, value, value_norm, value_cross_norm, dataset)] %>%
#     head()
# 
#   # Split datasets again
#   JYNR210602.sce2.norm.case <- ds1.sce4[dataset == "JYNR210602" & SampleId %in% rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,]),]
#   JYNR220601.sce2.norm.placebo <- ds1.sce4[dataset == "JYNR220601" & SampleId %in% rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]),]
# 
#   # Sanity check (should be 0)
#   length(unique(JYNR210602.sce2.norm.case$SampleId)) - dim(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 3,])[1] # should be 0
#   length(unique(JYNR220601.sce2.norm.placebo$SampleId)) - dim(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 1,])[1] # should be 0
# 
#   # Change to wide format and remove metadata
#   JYNR210602.sce2.norm.case <- JYNR210602.sce2.norm.case[order(SampleId), .(SampleId, variable, value_cross_norm)] %>%
#     dcast(SampleId ~ variable, value.var = "value_cross_norm")
#   JYNR210602.sce2.norm.case <- log2(as.data.frame(JYNR210602.sce2.norm.case)[,colnames(JYNR210602.sce2.norm.case) %in% colnames(log_JYNR210602)])
#   rownames(JYNR210602.sce2.norm.case) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,])
#   meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(JYNR210602.sce2.norm.case)]
# 
#   JYNR220601.sce2.norm.placebo <- JYNR220601.sce2.norm.placebo[order(SampleId), .(SampleId, variable, value_cross_norm)] %>%
#     dcast(SampleId ~ variable, value.var = "value_cross_norm")
#   JYNR220601.sce2.norm.placebo <- log2(as.data.frame(JYNR220601.sce2.norm.placebo)[,colnames(JYNR220601.sce2.norm.placebo) %in% colnames(log_JYNR220601)])
#   rownames(JYNR220601.sce2.norm.placebo) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,])
#   meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(JYNR220601.sce2.norm.placebo)]
  
  
  
## JYNR210602, SCENARIO 5
  # Derive reference values and merge
  ds1.sce5 <- rbind(ds1.long.nocont, ds2.long.nonash)
  ds1.sce5[SampleType %in% c('Sample', 'QC'), ref.sample.cross := median(value), .(variable)]
  
  # Apply normalization
  sample_norm(ds1.sce5, v = "value", ref = "ref.sample.cross", ext = 6)
  ds1.sce5[, value_cross_norm := value * sf6]
  ds1.sce5[, .(SampleId, variable, value, value_norm, value_cross_norm, dataset)] %>% 
    head()
  
  # Split datasets again
  JYNR210602.sce5.norm.case <- ds1.sce5[dataset == "JYNR210602" & SampleId %in% rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,]),]
  JYNR220601.sce5.norm.placebo <- ds1.sce5[dataset == "JYNR220601" & SampleId %in% rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]),]
  
  # Sanity check (should be 0)
  length(unique(JYNR210602.sce5.norm.case$SampleId)) - dim(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 3,])[1] # should be 0
  length(unique(JYNR220601.sce5.norm.placebo$SampleId)) - dim(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 1,])[1] # should be 0
  
  # Change to wide format and remove metadata
  JYNR210602.sce5.norm.case <- JYNR210602.sce5.norm.case[order(SampleId), .(SampleId, variable, value_cross_norm)] %>% 
    dcast(SampleId ~ variable, value.var = "value_cross_norm")
  JYNR210602.sce5.norm.case <- log2(as.data.frame(JYNR210602.sce5.norm.case)[,colnames(JYNR210602.sce5.norm.case) %in% colnames(log_JYNR210602)])
  rownames(JYNR210602.sce5.norm.case) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,])
  meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(JYNR210602.sce5.norm.case)]
  
  JYNR220601.sce5.norm.placebo <- JYNR220601.sce5.norm.placebo[order(SampleId), .(SampleId, variable, value_cross_norm)] %>% 
    dcast(SampleId ~ variable, value.var = "value_cross_norm")
  JYNR220601.sce5.norm.placebo <- log2(as.data.frame(JYNR220601.sce5.norm.placebo)[,colnames(JYNR220601.sce5.norm.placebo) %in% colnames(log_JYNR220601)])
  rownames(JYNR220601.sce5.norm.placebo) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,])
  meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(JYNR220601.sce5.norm.placebo)]
  
  
  
  
## JYNR220601, SCENARIO 1
  # Derive reference values and merge
  ds2.sce1 <- rbind(ds2.long.case, ds1.long.placebo)
  ds2.sce1[SampleType %in% c('Sample', 'QC'), ref.sample.cross := median(value), .(variable)]
  
  # Apply normalization
  sample_norm(ds2.sce1, v = "value", ref = "ref.sample.cross", ext = 6)
  ds2.sce1[, value_cross_norm := value * sf6]
  ds2.sce1[, .(SampleId, variable, value, value_norm, value_cross_norm, dataset)] %>% 
    head()
  
  # Split datasets again
  JYNR220601.sce1.norm.case <- ds2.sce1[dataset == "JYNR220601" & SampleId %in% rownames(log_JYNR220601),]
  JYNR210602.sce1.norm.placebo <- ds2.sce1[dataset == "JYNR210602" & SampleId %in% rownames(log_JYNR210602),]
  
  # Sanity check (should be 0)
  length(unique(JYNR220601.sce1.norm.case$SampleId)) - dim(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 2,])[1] # should be 0
  length(unique(JYNR210602.sce1.norm.placebo$SampleId)) - dim(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 1,])[1] # should be 0
  
  # Change to wide format and remove metadata
  JYNR220601.sce1.norm.case <- JYNR220601.sce1.norm.case[order(SampleId), .(SampleId, variable, value_cross_norm)] %>% 
    dcast(SampleId ~ variable, value.var = "value_cross_norm")
  JYNR220601.sce1.norm.case <- log2(as.data.frame(JYNR220601.sce1.norm.case)[,colnames(JYNR220601.sce1.norm.case) %in% colnames(log_JYNR220601)])
  rownames(JYNR220601.sce1.norm.case) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,])
  meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(JYNR220601.sce1.norm.case)]
  
  JYNR210602.sce1.norm.placebo <- JYNR210602.sce1.norm.placebo[order(SampleId), .(SampleId, variable, value_cross_norm)] %>% 
    dcast(SampleId ~ variable, value.var = "value_cross_norm")
  JYNR210602.sce1.norm.placebo <- log2(as.data.frame(JYNR210602.sce1.norm.placebo)[,colnames(JYNR210602.sce1.norm.placebo) %in% colnames(log_JYNR210602)])
  rownames(JYNR210602.sce1.norm.placebo) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,])
  meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(JYNR210602.sce1.norm.placebo)]
  
  
  
# ## JYNR220601, SCENARIO 2
#   # Derive reference values and merge
#   ds2.sce4 <- rbind(ds2.long.all.cases, ds1.long)
#   ds2.sce4[SampleType %in% c('Sample', 'QC'), ref.sample.cross := median(value), .(variable)]
# 
#   # Apply normalization
#   sample_norm(ds2.sce4, v = "value", ref = "ref.sample.cross", ext = 6)
#   ds2.sce4[, value_cross_norm := value * sf6]
#   ds2.sce4[, .(SampleId, variable, value, value_norm, value_cross_norm, dataset)] %>%
#     head()
# 
#   # Split datasets again
#   JYNR220601.sce2.norm.case <- ds2.sce4[dataset == "JYNR220601" & SampleId %in% rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,]),]
#   JYNR210602.sce2.norm.placebo <- ds2.sce4[dataset == "JYNR210602" & SampleId %in% rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,]),]
# 
#   # Sanity check (should be 0)
#   length(unique(JYNR220601.sce2.norm.case$SampleId)) - dim(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 2,])[1] # should be 0
#   length(unique(JYNR210602.sce2.norm.placebo$SampleId)) - dim(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 1,])[1] # should be 0
# 
#   # Change to wide format and remove metadata
#   JYNR220601.sce2.norm.case <- JYNR220601.sce2.norm.case[order(SampleId), .(SampleId, variable, value_cross_norm)] %>%
#     dcast(SampleId ~ variable, value.var = "value_cross_norm")
#   JYNR220601.sce2.norm.case <- log2(as.data.frame(JYNR220601.sce2.norm.case)[,colnames(JYNR220601.sce2.norm.case) %in% colnames(log_JYNR220601)])
#   rownames(JYNR220601.sce2.norm.case) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,])
#   meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(JYNR220601.sce2.norm.case)]
# 
#   JYNR210602.sce2.norm.placebo <- JYNR210602.sce2.norm.placebo[order(SampleId), .(SampleId, variable, value_cross_norm)] %>%
#     dcast(SampleId ~ variable, value.var = "value_cross_norm")
#   JYNR210602.sce2.norm.placebo <- log2(as.data.frame(JYNR210602.sce2.norm.placebo)[,colnames(JYNR210602.sce2.norm.placebo) %in% colnames(log_JYNR210602)])
#   rownames(JYNR210602.sce2.norm.placebo) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,])
#   meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(JYNR210602.sce2.norm.placebo)]
    
  
  
## JYNR220601, SCENARIO 5
  # Derive reference values and merge
  ds2.sce5 <- rbind(ds2.long.nocont, ds1.long.nonash)
  ds2.sce5[SampleType %in% c('Sample', 'QC'), ref.sample.cross := median(value), .(variable)]
  
  # Apply normalization
  sample_norm(ds2.sce5, v = "value", ref = "ref.sample.cross", ext = 6)
  ds2.sce5[, value_cross_norm := value * sf6]
  ds2.sce5[, .(SampleId, variable, value, value_norm, value_cross_norm, dataset)] %>% 
    head()
  
  # Split datasets again
  JYNR220601.sce5.norm.case <- ds2.sce5[dataset == "JYNR220601" & SampleId %in% rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,]),]
  JYNR210602.sce5.norm.placebo <- ds2.sce5[dataset == "JYNR210602" & SampleId %in% rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,]),]
  
  # Sanity check (should be 0)
  length(unique(JYNR220601.sce5.norm.case$SampleId)) - dim(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 2,])[1] # should be 0
  length(unique(JYNR210602.sce5.norm.placebo$SampleId)) - dim(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 1,])[1] # should be 0
  
  # Change to wide format and remove metadata
  JYNR220601.sce5.norm.case <- JYNR220601.sce5.norm.case[order(SampleId), .(SampleId, variable, value_cross_norm)] %>% 
    dcast(SampleId ~ variable, value.var = "value_cross_norm")
  JYNR220601.sce5.norm.case <- log2(as.data.frame(JYNR220601.sce5.norm.case)[,colnames(JYNR220601.sce5.norm.case) %in% colnames(log_JYNR220601)])
  rownames(JYNR220601.sce5.norm.case) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,])
  meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(JYNR220601.sce5.norm.case)]
  
  JYNR210602.sce5.norm.placebo <- JYNR210602.sce5.norm.placebo[order(SampleId), .(SampleId, variable, value_cross_norm)] %>% 
    dcast(SampleId ~ variable, value.var = "value_cross_norm")
  JYNR210602.sce5.norm.placebo <- log2(as.data.frame(JYNR210602.sce5.norm.placebo)[,colnames(JYNR210602.sce5.norm.placebo) %in% colnames(log_JYNR210602)])
  rownames(JYNR210602.sce5.norm.placebo) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,])
  meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(JYNR210602.sce5.norm.placebo)]
  


  
  
################################################################################
## Output dfs for batch effect correction (limma)
# Change to wide format and remove metadata
# JYNR210602, All groups, steps 1-4
ds1.wide <- ds1.long[SampleType == "Sample" & SampleId %in% rownames(log_JYNR210602)]
ds1.wide <- ds1.wide[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds1.wide <- log2(as.data.frame(ds1.wide)[,colnames(ds1.wide) %in% colnames(log_JYNR210602)])
rownames(ds1.wide) <- rownames(log_JYNR210602)
meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(ds1.wide)]

# JYNR210602, case, steps 1-4
ds1.wide.case <- ds1.long.case[SampleType == "Sample" & SampleId %in% rownames(log_JYNR210602)]
ds1.wide.case <- ds1.wide.case[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds1.wide.case <- log2(as.data.frame(ds1.wide.case)[,colnames(ds1.wide.case) %in% colnames(log_JYNR210602)])
rownames(ds1.wide.case) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,])
meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(ds1.wide.case)]

# JYNR210602, all cases, steps 1-4
ds1.wide.all.cases <- ds1.long.all.cases[SampleType == "Sample" & SampleId %in% rownames(log_JYNR210602)]
ds1.wide.all.cases <- ds1.wide.all.cases[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds1.wide.all.cases <- log2(as.data.frame(ds1.wide.all.cases)[,colnames(ds1.wide.all.cases) %in% colnames(log_JYNR210602)])
rownames(ds1.wide.all.cases) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group %in% c(3,4,5,6),])
meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(ds1.wide.all.cases)]

# JYNR210602, placebo, steps 1-4
ds1.wide.placebo <- ds1.long.placebo[SampleType == "Sample" & SampleId %in% rownames(log_JYNR210602)]
ds1.wide.placebo <- ds1.wide.placebo[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds1.wide.placebo <- log2(as.data.frame(ds1.wide.placebo)[,colnames(ds1.wide.placebo) %in% colnames(log_JYNR210602)])
rownames(ds1.wide.placebo) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,])
meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(ds1.wide.placebo)]

# JYNR210602, no NASH, steps 1-4
ds1.wide.nonash <- ds1.long.nonash[SampleType == "Sample" & SampleId %in% rownames(log_JYNR210602)]
ds1.wide.nonash <- ds1.wide.nonash[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds1.wide.nonash <- log2(as.data.frame(ds1.wide.nonash)[,colnames(ds1.wide.nonash) %in% colnames(log_JYNR210602)])
rownames(ds1.wide.nonash) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group %in% c(1,2,4,5,6),])
meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(ds1.wide.nonash)]

# JYNR210602, no control, steps 1-4
ds1.wide.nocont <- ds1.long.nocont[SampleType == "Sample" & SampleId %in% rownames(log_JYNR210602)]
ds1.wide.nocont <- ds1.wide.nocont[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds1.wide.nocont <- log2(as.data.frame(ds1.wide.nocont)[,colnames(ds1.wide.nocont) %in% colnames(log_JYNR210602)])
rownames(ds1.wide.nocont) <- rownames(log_JYNR210602[meta_mouse_JYNR210602$Group %in% c(2,3,4,5,6),])
meta_mouse_JYNR210602$Treatment[meta_mouse_JYNR210602$TubeID %in% rownames(ds1.wide.nocont)]




# JYNR220601, All groups, steps 1-4
ds2.wide <- ds2.long[SampleType == "Sample" & SampleId %in% rownames(log_JYNR220601)]
ds2.wide <- ds2.wide[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds2.wide <- log2(as.data.frame(ds2.wide)[,colnames(ds2.wide) %in% colnames(log_JYNR220601)])
rownames(ds2.wide) <- rownames(log_JYNR220601)
meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(ds2.wide)]

# JYNR220601, case, steps 1-4
ds2.wide.case <- ds2.long.case[SampleType == "Sample" & SampleId %in% rownames(log_JYNR220601)]
ds2.wide.case <- ds2.wide.case[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds2.wide.case <- log2(as.data.frame(ds2.wide.case)[,colnames(ds2.wide.case) %in% colnames(log_JYNR220601)])
rownames(ds2.wide.case) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,])
meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(ds2.wide.case)]

# JYNR220601, all cases, steps 1-4
ds2.wide.all.cases <- ds2.long.all.cases[SampleType == "Sample" & SampleId %in% rownames(log_JYNR220601)]
ds2.wide.all.cases <- ds2.wide.all.cases[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds2.wide.all.cases <- log2(as.data.frame(ds2.wide.all.cases)[,colnames(ds2.wide.all.cases) %in% colnames(log_JYNR220601)])
rownames(ds2.wide.all.cases) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group %in% c(2,3),])
meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(ds2.wide.all.cases)]

# JYNR220601, placebo, steps 1-4
ds2.wide.placebo <- ds2.long.placebo[SampleType == "Sample" & SampleId %in% rownames(log_JYNR220601)]
ds2.wide.placebo <- ds2.wide.placebo[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds2.wide.placebo <- log2(as.data.frame(ds2.wide.placebo)[,colnames(ds2.wide.placebo) %in% colnames(log_JYNR220601)])
rownames(ds2.wide.placebo) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,])
meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(ds2.wide.placebo)]

# JYNR220601, no NASH, steps 1-4
ds2.wide.nonash <- ds2.long.nonash[SampleType == "Sample" & SampleId %in% rownames(log_JYNR220601)]
ds2.wide.nonash <- ds2.wide.nonash[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds2.wide.nonash <- log2(as.data.frame(ds2.wide.nonash)[,colnames(ds2.wide.nonash) %in% colnames(log_JYNR220601)])
rownames(ds2.wide.nonash) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group %in% c(1,3),])
meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(ds2.wide.nonash)]

# JYNR220601, no control, steps 1-4
ds2.wide.nocont <- ds2.long.nocont[SampleType == "Sample" & SampleId %in% rownames(log_JYNR220601)]
ds2.wide.nocont <- ds2.wide.nocont[order(SampleId), .(SampleId, variable, value)] %>% 
  dcast(SampleId ~ variable, value.var = "value")
ds2.wide.nocont <- log2(as.data.frame(ds2.wide.nocont)[,colnames(ds2.wide.nocont) %in% colnames(log_JYNR220601)])
rownames(ds2.wide.nocont) <- rownames(log_JYNR220601[meta_mouse_JYNR220601$Group %in% c(2,3),])
meta_mouse_JYNR220601$Treatment[meta_mouse_JYNR220601$`ANIMAL#` %in% rownames(ds2.wide.nocont)]




################################################################################
## Cross-normalisation (batch effect correction) using limma package
library(limma)

# Merge datasets
ex1 <- t(log_JYNR210602)
ex2 <- t(log_JYNR220601)
ex <- merge(ex1, ex2, by = "row.names", all = TRUE)
rownames(ex) <- ex[,1]
ex <- ex[-1]
batch <- c(rep("JYNR210602",dim(ex1)[2]),rep("JYNR220601",dim(ex2)[2]))


# if (!is.null(batch)) {
#   batch <- as.factor(batch)
#   contrasts(batch) <- contr.sum(levels(batch))
#   batch <- model.matrix(~batch)[, -1, drop = FALSE]
# }
# x <- ex
# design <- matrix(1,ncol(x),1)
# X.batch <- cbind(batch)
# fit <- lmFit(x, cbind(design, X.batch))
# beta <- fit$coefficients[, -(1:ncol(design)), drop = FALSE]
# beta[is.na(beta)] <- 0
# as.matrix(x) - beta %*% t(X.batch)

# Remove batch effects
ex_rm <- removeBatchEffect(ex,batch)

# Split dfs
limma_JYNR210602 <- as.data.frame(t(ex_rm[,batch == "JYNR210602"]))
limma_JYNR220601 <- as.data.frame(t(ex_rm[,batch == "JYNR220601"]))



################################################################################
## Cross-normalisation (batch effect correction) using limma package
# For case:placebo pairings

## JYNR210602, SCENARIO 1, STEPS 1-4 = scenario 1
# Merge datasets
ex1 <- t(ds1.wide.case)
ex2 <- t(ds2.wide.placebo)
ex <- merge(ex1, ex2, by = "row.names", all = TRUE)
rownames(ex) <- ex[,1]
ex <- ex[-1]
batch <- c(rep("JYNR210602",dim(ex1)[2]),rep("JYNR220601",dim(ex2)[2]))

# Remove batch effects
ex_rm <- removeBatchEffect(ex,batch)

# Split dfs
JYNR210602.sce1.limma.case <- as.data.frame(t(ex_rm[,batch == "JYNR210602"]))
JYNR220601.sce1.limma.placebo <- as.data.frame(t(ex_rm[,batch == "JYNR220601"]))


# ## JYNR210602, SCENARIO 4, STEPS 1-4 = scenario 2
# # Merge datasets
# ex1 <- t(ds1.wide.all.cases)
# ex2 <- t(ds2.wide)
# ex <- merge(ex1, ex2, by = "row.names", all = TRUE)
# rownames(ex) <- ex[,1]
# ex <- ex[-1]
# batch <- c(rep("JYNR210602",dim(ex1)[2]),rep("JYNR220601",dim(ex2)[2]))
# 
# # Remove batch effects
# ex_rm <- removeBatchEffect(ex,batch)
# 
# # Split dfs
# JYNR210602.sce2.limma.case <- as.data.frame(t(ex_rm[,batch == "JYNR210602"])) %>% 
#   filter(rownames(ds1.wide.all.cases) %in% rownames(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 3,]))
# JYNR220601.sce2.limma.placebo <- as.data.frame(t(ex_rm[,batch == "JYNR220601"])) %>% 
#   filter(rownames(ds2.wide) %in% rownames(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 1,]))


## JYNR210602, SCENARIO 4, STEPS 1-4 = scenario 5
# Merge datasets
ex1 <- t(ds1.wide.nocont)
ex2 <- t(ds2.wide.nonash)
ex <- merge(ex1, ex2, by = "row.names", all = TRUE)
rownames(ex) <- ex[,1]
ex <- ex[-1]
batch <- c(rep("JYNR210602",dim(ex1)[2]),rep("JYNR220601",dim(ex2)[2]))

# Remove batch effects
ex_rm <- removeBatchEffect(ex,batch)

# Split dfs
JYNR210602.sce5.limma.case <- as.data.frame(t(ex_rm[,batch == "JYNR210602"])) %>% 
  filter(rownames(ds1.wide.nocont) %in% rownames(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 3,]))
JYNR220601.sce5.limma.placebo <- as.data.frame(t(ex_rm[,batch == "JYNR220601"])) %>% 
  filter(rownames(ds2.wide.nonash) %in% rownames(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 1,]))



## JYNR220601, SCENARIO 1, STEPS 1-4 = scenario 1
# Merge datasets
ex1 <- t(ds2.wide.case)
ex2 <- t(ds1.wide.placebo)
ex <- merge(ex1, ex2, by = "row.names", all = TRUE)
rownames(ex) <- ex[,1]
ex <- ex[-1]
batch <- c(rep("JYNR220601",dim(ex1)[2]),rep("JYNR210602",dim(ex2)[2]))

# Remove batch effects
ex_rm <- removeBatchEffect(ex,batch)

# Split dfs
JYNR220601.sce1.limma.case <- as.data.frame(t(ex_rm[,batch == "JYNR220601"]))
JYNR210602.sce1.limma.placebo <- as.data.frame(t(ex_rm[,batch == "JYNR210602"]))


# ## JYNR220601, SCENARIO 4, STEPS 1-4 = scenario 2
# # Merge datasets
# ex1 <- t(ds2.wide.all.cases)
# ex2 <- t(ds1.wide)
# ex <- merge(ex1, ex2, by = "row.names", all = TRUE)
# rownames(ex) <- ex[,1]
# ex <- ex[-1]
# batch <- c(rep("JYNR220601",dim(ex1)[2]),rep("JYNR210602",dim(ex2)[2]))
# 
# # Remove batch effects
# ex_rm <- removeBatchEffect(ex,batch)
# 
# # Split dfs
# JYNR220601.sce2.limma.case <- as.data.frame(t(ex_rm[,batch == "JYNR220601"])) %>% 
#   filter(rownames(ds2.wide.all.cases) %in% rownames(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 2,]))
# JYNR210602.sce2.limma.placebo <- as.data.frame(t(ex_rm[,batch == "JYNR210602"])) %>% 
#   filter(rownames(ds1.wide) %in% rownames(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 1,]))


## JYNR220601, SCENARIO 4, STEPS 1-4 = scenario 5
# Merge datasets
ex1 <- t(ds2.wide.nocont)
ex2 <- t(ds1.wide.nonash)
ex <- merge(ex1, ex2, by = "row.names", all = TRUE)
rownames(ex) <- ex[,1]
ex <- ex[-1]
batch <- c(rep("JYNR220601",dim(ex1)[2]),rep("JYNR210602",dim(ex2)[2]))

# Remove batch effects
ex_rm <- removeBatchEffect(ex,batch)

# Split dfs
JYNR220601.sce5.limma.case <- as.data.frame(t(ex_rm[,batch == "JYNR220601"])) %>% 
  filter(rownames(ds2.wide.nocont) %in% rownames(meta_mouse_JYNR220601[meta_mouse_JYNR220601$Group == 2,]))
JYNR210602.sce5.limma.placebo <- as.data.frame(t(ex_rm[,batch == "JYNR210602"])) %>% 
  filter(rownames(ds1.wide.nonash) %in% rownames(meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 1,]))




################################################################################
## Compile all case/placebo paired normalisation dfs into one list
case_control_norm_dfs <- list(list(JYNR210602.sce1.norm.case,
                                   JYNR220601.sce1.norm.placebo),
                              list(JYNR220601.sce1.norm.case,
                                   JYNR210602.sce1.norm.placebo),
                              list(JYNR210602.sce1.limma.case,
                                   JYNR220601.sce1.limma.placebo),
                              list(JYNR220601.sce1.limma.case,
                                   JYNR210602.sce1.limma.placebo),
                              list(JYNR210602.sce5.norm.case,
                                   JYNR220601.sce5.norm.placebo),
                              list(JYNR220601.sce5.norm.case,
                                   JYNR210602.sce5.norm.placebo),
                              list(JYNR210602.sce5.limma.case,
                                   JYNR220601.sce5.limma.placebo),
                              list(JYNR220601.sce5.limma.case,
                                   JYNR210602.sce5.limma.placebo))



################################################################################
## Calibrator median normalisation
# If I had a lot of calibrator samples in each study,
# it might make the most sense to compute scale factors for each SOMAmer.
# However, since I only have 3 calibrator samples in each study,
# I am inclined to believe that differences on a SOMAmer-basis
# are the result of natural sample variance, not systematic differences.
# Thus, I do not believe they reflect the actual inter-study differences.
# Hence, this form of normalisation would not reduce inter-study noise.
# We therefore have two options - using a median for all SOMAmers,
# or computing a median within each dilution bin.
# (Let me just quickly assess the difference in absolute values between the bins)
get_meds <- function(df){
  low <- df[,colnames(df) %in% SomaScanAnnotation$Seq[SomaScanAnnotation$Dilution == "20%"]]
  med <- df[,colnames(df) %in% SomaScanAnnotation$Seq[SomaScanAnnotation$Dilution == "0.5%"]]
  high <- df[,colnames(df) %in% SomaScanAnnotation$Seq[SomaScanAnnotation$Dilution == "0.005%"]]
  low <- sapply(1:ncol(low), function(x) median(low[,x]))
  med <- sapply(1:ncol(med), function(x) median(med[,x]))
  high <- sapply(1:ncol(high), function(x) median(high[,x]))
  # xlim <- c(min(c(low,med,high)),max(c(low,med,high)))
  # layout(mat = matrix(c(1,2,3), nrow = 3, ncol = 1, byrow=TRUE))
  # par(oma=c(0,0,3,1))
  # par(mar=c(4.5,5.5,2.25,0))
  # hist(low, breaks=50, xlab = "Median RFU per aptamer", xlim = xlim, main = "Dilution = 20%")
  # abline(v = median(low),lwd=2.5,col="red")
  # hist(med, breaks=50, xlab = "Median RFU per aptamer", xlim = xlim, main = "Dilution = 0.5%")
  # abline(v = median(med),lwd=2.5,col="red")
  # hist(high, breaks=50, xlab = "Median RFU per aptamer", xlim = xlim, main = "Dilution = 0.005%")
  # abline(v = median(high),lwd=2.5,col="red")
  
  meds <- list(low, med, high)
  return(meds)
}
# *tested some dfs*
# yeah, let's work with scale factors for individual dilution bins
# they are on different scales
# and let's work with log2-scaled data


meds_1 <- get_meds(log2(cal_JYNR210602))
meds_2 <- get_meds(log2(cal_JYNR220601))

# For all aptamers (not differentiating btw. dilution bins)
med_1 <- unlist(meds_1); med_2 <- unlist(meds_2)
ratios <- sapply(1:length(med_1), function(x) med_1[x]/med_2[x])
scale_factor <- median(ratios)
# This yields a scale factor of 1.0002
# One SF
cal_one_JYNR210602 <- log_JYNR210602 * scale_factor

# For the three dilution bins
scale_factors <- c()
for (dil in 1:3){
  med_1 <- meds_1[[dil]]; med_2 <- meds_2[[dil]]
  ratios <- sapply(1:length(med_1), function(x) med_1[x]/med_2[x])
  scale_factor <- median(ratios)
  scale_factors <- c(scale_factors, scale_factor)
}
# This yields the scale factors 1.001, 0.996, and 0.996
# These will not make a difference.
cal_three_JYNR210602 <- log_JYNR210602
cal_three_JYNR210602[colnames(log_JYNR210602) %in% SomaScanAnnotation$Seq[SomaScanAnnotation$Dilution %in% c("20%")]] <- log_JYNR210602[,colnames(log_JYNR210602) %in% SomaScanAnnotation$Seq[SomaScanAnnotation$Dilution %in% c("20%")]] * scale_factors[1]
cal_three_JYNR210602[colnames(log_JYNR210602) %in% SomaScanAnnotation$Seq[SomaScanAnnotation$Dilution %in% c("0.5%")]] <- log_JYNR210602[,colnames(log_JYNR210602) %in% SomaScanAnnotation$Seq[SomaScanAnnotation$Dilution %in% c("0.5%")]] * scale_factors[2]
cal_three_JYNR210602[colnames(log_JYNR210602) %in% SomaScanAnnotation$Seq[SomaScanAnnotation$Dilution %in% c("0.005")]] <- log_JYNR210602[,colnames(log_JYNR210602) %in% SomaScanAnnotation$Seq[SomaScanAnnotation$Dilution %in% c("0.005")]] * scale_factors[3]

# Per aptamer-basis
med_1 <- unlist(meds_1); med_2 <- unlist(meds_2)
ratios <- sapply(1:length(med_1), function(x) med_1[x]/med_2[x])
min(ratios); max(ratios)
# Between 0.55 and 1.37
# hist(ratios)
sum(ratios < 0.9) # 19
sum(ratios > 1.1) # 22
cal_indi_JYNR210602 <- log_JYNR210602 * ratios
# This will likely not make a difference either.
# hist(ratios, main = "7212 scale factors", xlab = "Scale factors")

# The conclusion seems to be that the calibrator samples are quite similar.
# They have, after all, been calibrated to the same platescale references.
# It seems like the cross-normalisation on a per aptamer-basis might be feasible.


################################################################################
## Remove unnecessary items from the environment
print("Cross-normalisation complete.")
print("Removing unnecessary items from the environment...")
rm(list = c("unnorm1", "unnorm2", "hybNorm.medNormInt.plateScale_orig", "crossnorm",
            "ds.merged",
            "ds1", "ds1.sce1", "ds1.sce5",
            "ds1.long", "ds1.long.case", "ds1.long.all.cases", "ds1.long.placebo",
            "ds1.wide", "ds1.wide.case", "ds1.wide.all.cases", "ds1.wide.placebo",
            "ds1.long.nocont", "ds1.long.nonash", "ds1.wide.nocont", "ds1.wide.nonash",
            "ds2", "ds2.sce1", "ds2.sce5",
            "ds2.long", "ds2.long.case", "ds2.long.all.cases", "ds2.long.placebo",
            "ds2.wide", "ds2.wide.case", "ds2.wide.all.cases", "ds2.wide.placebo",
            "ds2.long.nocont", "ds2.long.nonash", "ds2.wide.nocont", "ds2.wide.nonash",
            "ex",
            "ex_rm",
            "ex1",
            "ex2",
            "batch",
            "sample_norm",
            "meta",
            "seq_cols",
            "JYNR210602.sce1.norm.case", "JYNR210602.sce1.norm.placebo",
            "JYNR220601.sce1.norm.case", "JYNR220601.sce1.norm.placebo",
            "JYNR210602.sce1.limma.case", "JYNR210602.sce1.limma.placebo",
            "JYNR220601.sce1.limma.case", "JYNR220601.sce1.limma.placebo",
            # "JYNR210602.sce2.norm.case", "JYNR210602.sce2.norm.placebo",
            # "JYNR220601.sce2.norm.case", "JYNR220601.sce2.norm.placebo",
            # "JYNR210602.sce2.limma.case", "JYNR210602.sce2.limma.placebo",
            # "JYNR220601.sce2.limma.case", "JYNR220601.sce2.limma.placebo",
            # "JYNR210602.sce3.log.case", "JYNR210602.sce3.log.placebo",
            # "JYNR220601.sce3.log.case", "JYNR220601.sce3.log.placebo",
            # "JYNR210602.sce3.limma.case", "JYNR210602.sce3.limma.placebo",
            # "JYNR220601.sce3.limma.case", "JYNR220601.sce3.limma.placebo",
            # "JYNR210602.sce4.log.case", "JYNR210602.sce4.log.placebo",
            # "JYNR220601.sce4.log.case", "JYNR220601.sce4.log.placebo",
            # "JYNR210602.sce4.limma.case", "JYNR210602.sce4.limma.placebo",
            # "JYNR220601.sce4.limma.case", "JYNR220601.sce4.limma.placebo",
            "JYNR210602.sce5.norm.case", "JYNR210602.sce5.norm.placebo",
            "JYNR220601.sce5.norm.case", "JYNR220601.sce5.norm.placebo",
            "JYNR210602.sce5.limma.case", "JYNR210602.sce5.limma.placebo",
            "JYNR220601.sce5.limma.case", "JYNR220601.sce5.limma.placebo",
            "meds_1", "meds_2", "med_1", "med_2",
            "dil", "i", "get_meds",
            "ratios", "scale_factor", "scale_factors"
            ))

