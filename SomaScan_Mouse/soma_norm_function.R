## Initialize document
library(tidyverse); library(magrittr); library(data.table); library(SomaDataIO)
#setwd("~/JDJG/clean")


## Load test data
# unnorm1 <- read_adat("~/NNEDL/masterprojectxjdjg/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.adat")
# unnorm2 <- read_adat("~/NNEDL/masterprojectxjdjg/raw/Mouse_studies/JYNR220601_GUS2022_898/SS-2231349_v4.1_other.adat")
# hybNorm.medNormInt.plateScale_orig <- read_adat("~/NNEDL/masterprojectxjdjg/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.hybNorm.medNormInt.plateScale.adat")
# hybNorm.medNormInt.plateScale.medNormSMP_orig <- read_adat("~/NNEDL/masterprojectxjdjg/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.hybNorm.medNormInt.plateScale.medNormSMP.adat")
# 
# load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")

# Define function values for manual walkthrough
# adat <- unnorm
# annotation_df <- hybNorm.medNormInt.plateScale
# medNormSMP = FALSE


## Define the normalisation function
soma_norm <- function(adat, annotation_df, hybNorm = TRUE, medNormInt = TRUE,
                      plateScale = TRUE, medNormSMP = TRUE){
  # Note: annotation_df needs to be a normalised adat for the PlateScale_References
  
  
  ## Separate assay info from RFU values
  seq_cols <- names(adat) %>%  { .[grepl("^seq", .)] }
  assay_data <- adat[,!names(adat) %in% seq_cols]
  adat <- adat[,seq_cols]
  annotation <- getAnalyteInfo(annotation_df)
  
  
  ## Hybridisation normalisation
  if (hybNorm == TRUE){
    # Get column names of hybridization sequences with getAnalyteInfo function
    # (they will have Type == "Hybridization Control Elution")
    hyb_cols <- annotation$AptName[annotation$Type == "Hybridization Control Elution"]
    # Take median of each sequence across all samples to be the reference for each hyb sequence 
    references <- sapply(hyb_cols, function(x) median(adat[,x]))
    # For each sample, take a ratio of  hyb seq reference divided by the seq in the given sample
    # (there are twelve sequences, so should be twelve ratios)
    ratios <- matrix(nrow = dim(adat)[1], ncol = length(hyb_cols))
    for (i in 1:dim(ratios)[1]){
      for (j in 1:dim(ratios)[2])
        ratios[i,j] <- references[j] / adat[i,hyb_cols[j]]
    }
    # Take the median of ratios, this is the hyb scale factor
    scale_factors <- sapply(1:dim(ratios)[1], function(x) median(ratios[x,]))
    adat <- adat * scale_factors
    adat <- adat %>% round(1)
  }
  
  
  
  ## Intraplate median normalisation (just calibrators and buffers)
  if (medNormInt == TRUE){
    # Calibrator samples
    # Get the column names of each aptamer in the 20%, 0.5%, and 0.005% dilution groups
    dil_cols <- list(annotation$AptName[annotation$Dilution == "20"],
                     annotation$AptName[annotation$Dilution == "0.5"],
                     annotation$AptName[annotation$Dilution == "0.005"])
    cal_samples <- rownames(assay_data)[assay_data$SampleType == "Calibrator"]
    # Take the median of each analyte across calibrators on a given plate to be the intraplate calibrator reference
    references <- sapply(seq_cols, function(x) median(adat[cal_samples,x]))
    # For each sample, take a ratio of the intraplate calibrator reference divided 
    # by the value in calibrator replicate (a ratio for each analyte)
    ratios <- matrix(nrow = length(cal_samples), ncol = length(seq_cols))
    colnames(ratios) <- colnames(adat)
    for (i in 1:dim(ratios)[1]){
      for (j in 1:dim(ratios)[2]){
        ratios[i,j] <- references[j] / adat[cal_samples[i],j]
      }
    }
    for (i in 1:3){
      # Take the median of ratios of aptamers in a given dilution group, this is the intraplate scale factor
      scale_factors <- sapply(1:dim(ratios)[1], function(x) median(ratios[x,dil_cols[[i]]]))
      # Multiply this factor by the analytes in that calibrator in that dilution group
      adat[cal_samples,dil_cols[[i]]] <- adat[cal_samples,dil_cols[[i]]] * scale_factors
    }
    # Repeat for buffer samples
    # Get the column names of each aptamer in the 20%, 0.5%, and 0.005% dilution groups
    buf_samples <- rownames(assay_data)[assay_data$SampleType == "Buffer"]
    # Take the median of each analyte across calibrators on a given plate to be the intraplate calibrator reference
    references <- sapply(seq_cols, function(x) median(adat[buf_samples,x]))
    # For each sample, take a ratio of the intraplate calibrator reference divided 
    # by the value in calibrator replicate (a ratio for each analyte)
    ratios <- matrix(nrow = length(buf_samples), ncol = length(seq_cols))
    colnames(ratios) <- colnames(adat)
    for (i in 1:dim(ratios)[1]){
      for (j in 1:dim(ratios)[2])
        ratios[i,j] <- references[j] / adat[buf_samples[i],j]
    }
    for (i in 1:3){
      # Take the median of ratios of aptamers in a given dilution group, this is the intraplate scale factor
      scale_factors <- sapply(1:dim(ratios)[1], function(x) median(ratios[x,dil_cols[[i]]]))
      # Multiple this factor by the analytes in that calibrator in that dilution group
      adat[buf_samples,dil_cols[[i]]] <- adat[buf_samples,dil_cols[[i]]] * scale_factors
    }
    adat <- adat %>% round(1)
  }
  
  
  
  ## Plate scaling
  if (plateScale == TRUE){
    # Note: Only designed to work with data from a single plate, not multiple plates
    # Get the plate scale reference for each analyte
    # (this can be found with getAnalyteInfo in the `PlateScale_Reference` column 
    #  of any normalized data file of that sample type)
    references <- annotation$PlateScale_Reference # common reference for all matrices
    # Take the median of each analyte the now-normalized calibrators on a given plate
    cal_samples <- rownames(assay_data)[assay_data$SampleType == "Calibrator"]
    medians <- sapply(seq_cols, function(x) median(adat[cal_samples,x]))
    # Take a ratio of the plate scale reference divided by median of calibrators on the plate
    # (do this for each analyte)
    ratios <- sapply(1:length(seq_cols), function(x) references[x] / medians[x])
    # Take the median of ratios, this is the plate scale factor
    # (Round plate scale factor to seven decimal places)
    # platescale_factor <- round(median(ratios),7)
    platescale_factor <- median(ratios)
    # Multiply this factor by all samples/proteins on the plate
    adat <- round(adat * platescale_factor,1)
  }
  
  
  
  ## Interplate median normalisation (just samples)
  if (medNormSMP == TRUE){
    # perform Interplate median normalization
    # Get the column names of each aptamer in the 20%, 0.5%, and 0.005% dilution groups
    dil_cols <- list(annotation$AptName[annotation$Dilution == "20"],
                     annotation$AptName[annotation$Dilution == "0.5"],
                     annotation$AptName[annotation$Dilution == "0.005"])
    # Take the median of each analyte across all samples and across all plates. This is the sample median reference
    samples <- rownames(assay_data)[assay_data$SampleType == "Sample"]
    references <- sapply(seq_cols, function(x) median(adat[samples,x]))
    # For each sample, take a ratio of the sample median reference divided
    # by the value in the sample (do this for each analyte)
    ratios <- matrix(nrow = length(samples), ncol = length(references))
    colnames(ratios) <- colnames(adat[samples,])
    for (i in 1:dim(ratios)[1]){
      for (j in 1:dim(ratios)[2])
        ratios[i,j] <- references[j] / adat[samples[i],j]
    }
    # Take the median of ratios of aptamers in a given dilution group, this is the interplate median scale factor
    for (i in 1:3){
      # Take the median of ratios of aptamers in a given dilution group, this is the intraplate scale factor
      scale_factors <- sapply(1:dim(ratios)[1], function(x) median(ratios[x,dil_cols[[i]]]))
      # Multiple this factor by the analytes in that calibrator in that dilution group
      adat[samples,dil_cols[[i]]] <- adat[samples,dil_cols[[i]]] * scale_factors
    }
    adat <- adat %>% round(1)
    
  }
  
  adat <- merge(assay_data, adat, by = "row.names") # keep assay data, just in case
  rownames(adat) <- adat$Row.names; adat <- adat[,-1]
  
  return(adat)
}


################################################################################
## Cross-normalisation function
# As per the recommendations of SomaLogic
## Define the normalisation function
soma_crossnorm <- function(adat1, adat2, SomaScanAnnotation = SomaScanAnnotation){
  adat1 <- log_JYNR210602; adat2 <- log_JYNR220601
  merged_adat <- full_join(as.data.frame(adat1), as.data.frame(adat2))
  rownames(merged_adat) <- c(rownames(adat1), rownames(adat2))
  
  adat <- hybNorm.medNormInt.plateScale_orig
  seq_cols <- names(adat) %>%  { .[grepl("^seq", .)] }
  assay_data <- adat[,!names(adat) %in% seq_cols]
  adat <- adat[,seq_cols]
  annotation <- getAnalyteInfo(annotation_df)
  merged_adat <- adat
  
  # join the two dfs and simply run one norm step (identical to the one in the function above)
  # split the dfs and return a list

  ## Interplate median normalisation (just samples)
  # Get the column names of each aptamer in the 20%, 0.5%, and 0.005% dilution groups
  dil_cols <- list(SomaScanAnnotation$SeqId[SomaScanAnnotation$Dilution == "20%" & SomaScanAnnotation$SeqId %in% colnames(merged_adat)],
                   SomaScanAnnotation$SeqId[SomaScanAnnotation$Dilution == "0.5%" & SomaScanAnnotation$SeqId %in% colnames(merged_adat)],
                   SomaScanAnnotation$SeqId[SomaScanAnnotation$Dilution == "0.005%" & SomaScanAnnotation$SeqId %in% colnames(merged_adat)])
  # Take the median of each analyte across all samples and across all plates. This is the sample median reference
  references <- sapply(colnames(merged_adat), function(x) median(merged_adat[,x]))
  # For each sample, take a ratio of the sample median reference divided
  # by the value in the sample (do this for each analyte)
  ratios <- matrix(nrow = dim(merged_adat)[1], ncol = length(references))
  colnames(ratios) <- colnames(merged_adat)
  for (i in 1:dim(ratios)[1]){
    for (j in 1:dim(ratios)[2])
      ratios[i,j] <- references[j] / merged_adat[i,j]
  }
  # Take the median of ratios of aptamers in a given dilution group, this is the interplate median scale factor
  for (i in 1:3){
    # Take the median of ratios of aptamers in a given dilution group, this is the intraplate scale factor
    scale_factors <- sapply(1:dim(ratios)[1], function(x) median(ratios[x,dil_cols[[i]]]))
    # Multiple this factor by the analytes in that calibrator in that dilution group
    merged_adat[,dil_cols[[i]]] <- merged_adat[,dil_cols[[i]]] * scale_factors
  }
  merged_adat <- merged_adat %>% round(1)
  
  norm_adat1 <- merged_adat[rownames(merged_adat) %in% rownames(adat1),]
  norm_adat2 <- merged_adat[rownames(merged_adat) %in% rownames(adat2),]
  
  return (list(norm_adat1, norm_adat2))
 
}


soma_crossnorm <- function(adat1, adat2, annotation_df){
  # adat1 <- log_JYNR210602; adat2 <- log_JYNR220601
  norm_1 <- soma_norm(adat = adat1, annotation_df = annotation_df, medNormSMP = FALSE)
  norm_2 <- soma_norm(adat = adat2, annotation_df = annotation_df, medNormSMP = FALSE)
  
  # Join adat dataframes
  merged_adat <- full_join(as.data.frame(adat1), as.data.frame(adat2))
  rownames(merged_adat) <- c(rownames(adat1), rownames(adat2)) # keep row names for separation later
  
  # adat <- hybNorm.medNormInt.plateScale_orig
  seq_cols <- names(merged_adat) %>%  { .[grepl("^seq", .)] }
  assay_data <- merged_adat[,!names(merged_adat) %in% seq_cols]
  hybNorm.medNormInt.plateScale <- merged_adat[,seq_cols]
  annotation <- getAnalyteInfo(annotation_df)
  
  # join the two dfs and simply run one norm step (identical to the one in the function above)
  # split the dfs and return a list
  
  
  # perform Interplate median normalization
  # Get the column names of each aptamer in the 20%, 0.5%, and 0.005% dilution groups
  dil_cols <- list(annotation$AptName[annotation$Dilution == "20"],
                   annotation$AptName[annotation$Dilution == "0.5"],
                   annotation$AptName[annotation$Dilution == "0.005"])
  # Take the median of each analyte across all samples and across all plates. This is the sample median reference
  samples <- rownames(assay_data)[assay_data$SampleType == "Sample"]
  references <- sapply(seq_cols, function(x) median(hybNorm.medNormInt.plateScale[samples,x]))
  # For each sample, take a ratio of the sample median reference divided
  # by the value in the sample (do this for each analyte)
  ratios <- matrix(nrow = length(samples), ncol = length(references))
  colnames(ratios) <- colnames(hybNorm.medNormInt.plateScale[samples,])
  for (i in 1:dim(ratios)[1]){
    for (j in 1:dim(ratios)[2])
      ratios[i,j] <- references[j] / hybNorm.medNormInt.plateScale[samples[i],j]
  }
  # Take the median of ratios of aptamers in a given dilution group, this is the interplate median scale factor
  hybNorm.medNormInt.plateScale.medNormSMP <- hybNorm.medNormInt.plateScale
  for (i in 1:3){
    # Take the median of ratios of aptamers in a given dilution group, this is the intraplate scale factor
    scale_factors <- sapply(1:dim(ratios)[1], function(x) median(ratios[x,dil_cols[[i]]]))
    # Multiple this factor by the analytes in that calibrator in that dilution group
    hybNorm.medNormInt.plateScale.medNormSMP[samples,dil_cols[[i]]] <- hybNorm.medNormInt.plateScale.medNormSMP[samples,dil_cols[[i]]] * scale_factors
  }
  hybNorm.medNormInt.plateScale.medNormSMP <- hybNorm.medNormInt.plateScale.medNormSMP %>% round(1)

  crossnorm_1 <- hybNorm.medNormInt.plateScale.medNormSMP[rownames(hybNorm.medNormInt.plateScale.medNormSMP) %in% rownames(adat1),]
  assay_data_1 <- assay_data[rownames(assay_data) %in% rownames(adat1),]
  crossnorm_1 <- merge(assay_data_1, crossnorm_1, by = "row.names") # keep assay data, just in case
  rownames(crossnorm_1) <- crossnorm_1$Row.names; crossnorm_1 <- crossnorm_1[,-1]

  crossnorm_2 <- hybNorm.medNormInt.plateScale.medNormSMP[rownames(hybNorm.medNormInt.plateScale.medNormSMP) %in% rownames(adat2),]
  assay_data_2 <- assay_data[rownames(assay_data) %in% rownames(adat2),]
  crossnorm_2 <- merge(assay_data_2, crossnorm_2, by = "row.names") # keep assay data, just in case
  rownames(crossnorm_2) <- crossnorm_2$Row.names; crossnorm_2 <- crossnorm_2[,-1]
  
  return (list(crossnorm_1, crossnorm_2))
  
}



# ## Testing the function
# homemade_norm <- soma_norm(adat = unnorm,
#                   annotation_df = hybNorm.medNormInt.plateScale.medNormSMP_orig,
#                   medNormSMP = TRUE)
# fake_crossnorm <- soma_crossnorm(adat1 = unnorm1, adat2 = unnorm1, annotation_df = hybNorm.medNormInt.plateScale_orig)[[1]]
# 
# 
# 
# # Make adat comparable to fully wrangled data:
# adat <- crossnorm_1
# test <- adat[adat$SampleType == "Sample",]
# rownames(test) <- test$SampleId
# test <- test[rownames(test) %in% rownames(log_JYNR210602),
#              colnames(test) %in% colnames(log_JYNR210602)]
# test <- test[order(rownames(test)),]
# rownames(test) == rownames(log_JYNR210602)
# test <- log2(as.data.frame(test))
# 
# 
# # norm_adat1 = SomaLogic's 1-4 + mit nye 5'er script
# # Test mod din gamle 5'er
# #log_1_5_JYNR210602 <- log_JYNR210602
# #crossnorm_1_5_JYNR210602 <- test
# #crossnorm_1_5_JYNR210602_test2 <- test
# #norm_1_5_JYNR210602 <- test
# 
# # Checking equivalency:
# # soma_norm <- as.data.frame(hybNorm.medNormInt.plateScale.medNormSMP_orig[,seq_cols])
# # soma_norm <- as.data.frame(merged_adat)
# 
# test[1:5,1:5] == as.data.frame(log_JYNR210602)[1:5,1:5]
# test[1:5,1:5] == hybNorm.medNormInt.plateScale.medNormSMP_orig[1:5,1:5]
# as.data.frame(log_JYNR210602)[1:5,1:5]
# 
# crossnorm_1_5_JYNR210602[1:5,1:5]
# crossnorm_1_5_JYNR210602_test2[1:5,1:5]
# norm_1_5_JYNR210602[1:5,1:5]
# log_1_5_JYNR210602[1:5,1:5]
# 
# rownames(meta_mouse_JYNR210602)
# meta_mouse_JYNR210602$
# 
# sum(!homemade_norm == soma_norm) # only 14 values do not match (presumably due to rounding differences)