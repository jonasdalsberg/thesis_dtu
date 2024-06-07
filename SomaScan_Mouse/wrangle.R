################################################################################
library(SomaDataIO)

## Set working directory
setwd("~/JDJG/clean/")

### Mouse data
## Load data
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData") # load curated data (to avoid having to run everything each time)
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse.RData")

# Get annotation file
SomaScanAnnotation <- read.csv(paste0("~/NNEDL/somalogicannotation/curated/v4_1_7k/current/", list.files("~/NNEDL/somalogicannotation/curated/v4_1_7k/current/")))
SomaScanAnnotation$SeqId <- paste0("seq.",sub("\\-.*", "", SomaScanAnnotation$SeqId),".",sub('.+-(.+)', '\\1', SomaScanAnnotation$SeqId))

## Wrangle
# JYNR210602 (GUS2021-390-NN)
meta_mouse_JYNR210602 <- Phenotype_Mouse_JYNR210602 %>% 
  filter(TubeID %in% Phenotype_Mouse_JYNR210602_w_soma$Unique_ID)
meta_mouse_JYNR210602$id <- meta_mouse_JYNR210602$TubeID
meta_mouse_JYNR210602 <- meta_mouse_JYNR210602[order(meta_mouse_JYNR210602$id),]
for (i in c(16,31,58)){
  meta_mouse_JYNR210602[,i] <- as.character(meta_mouse_JYNR210602[,i])
}
rownames(meta_mouse_JYNR210602) <- meta_mouse_JYNR210602$id
meta_mouse_JYNR210602$dataset <- "SomaMouse1"
mouse_JYNR210602 <- Phenotype_Mouse_JYNR210602_w_soma %>% 
  select(-Unique_ID)
mouse_JYNR210602 <- mouse_JYNR210602[order(rownames(mouse_JYNR210602)),]
mouse_JYNR210602 <- mouse_JYNR210602[,-c(7225:7305)]
mouse_JYNR210602 <- mouse_JYNR210602[,!colnames(mouse_JYNR210602) %in% SomaScanAnnotation$SeqId[SomaScanAnnotation$Type == "Hybridization Control Elution"]]
# for (i in 1:dim(mouse_JYNR210602)[2]){
#   if (!(str_detect(colnames(mouse_JYNR210602)[i], "seq"))){
#     print(i)
#     print(colnames(mouse_JYNR210602)[i])
#   }
# }


# JYNR220601 (GUS2022-898-NN)
meta_mouse_JYNR220601 <- as.data.frame(Phenotype_Mouse_JYNR220601) %>% 
  filter(`ANIMAL#` %in% Phenotype_Mouse_JYNR220601_w_soma$Unique_ID) %>% 
  select(-Team)
meta_mouse_JYNR220601$id <- meta_mouse_JYNR220601$`ANIMAL#`
meta_mouse_JYNR220601 <- meta_mouse_JYNR220601[order(meta_mouse_JYNR220601$id),]
meta_mouse_JYNR220601$Group[meta_mouse_JYNR220601$Group == 8] <- 3
meta_mouse_JYNR220601$Diet[meta_mouse_JYNR220601$Diet == "D09100310-Ssniff (40% GAN)"] <- "D09100310 (GAN 40%)"
for (i in c(6,9:11,13:17,20:36,54:74)){
  meta_mouse_JYNR220601[,i] <- as.numeric(meta_mouse_JYNR220601[,i])
}
rownames(meta_mouse_JYNR220601) <- meta_mouse_JYNR220601$id
meta_mouse_JYNR220601$dataset <- "SomaMouse2"
mouse_JYNR220601 <- Phenotype_Mouse_JYNR220601_w_soma %>% 
  select(-Unique_ID)
mouse_JYNR220601 <- mouse_JYNR220601[order(rownames(mouse_JYNR220601)),]
mouse_JYNR220601 <- mouse_JYNR220601[,-c(7242:7316)]
# for (i in 1:dim(mouse_JYNR220601)[2]){
#   if (!(str_detect(colnames(mouse_JYNR220601)[i], "seq"))){
#     print(i)
#     print(colnames(mouse_JYNR220601)[i])
#   }
# }


# Align dataframes (and censor names)
meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 2,]$Treatment <- "Chow Treatment A"
meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 4,]$Treatment <- "DIO-NASH Treatment A"
meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 5,]$Treatment <- "DIO-NASH Treatment B"
meta_mouse_JYNR210602[meta_mouse_JYNR210602$Group == 6,]$Treatment <- "DIO-NASH Treatment C"
meta_mouse_JYNR220601 <- meta_mouse_JYNR220601 %>% 
  mutate(Model = case_when(
    Diet == "1324 (Altromin chow)" ~ "Healthy",
    Diet == "D09100310 (GAN 40%)" ~ "NASH"))
meta_mouse_JYNR220601 <- meta_mouse_JYNR220601 %>% 
  mutate(Treatment = case_when(
    Treatment == "Vehicle" ~ "DIO-NASH Vehicle",
    Treatment == "Semaglutide" ~ "DIO-NASH Treatment D",
    Treatment == "Lean chow vehicle" ~ "Chow Vehicle"))
meta_mouse_JYNR210602$Treatment_comp <- meta_mouse_JYNR210602$Treatment
meta_mouse_JYNR220601$Treatment_comp <- meta_mouse_JYNR220601$Treatment
meta_mouse_JYNR210602$Treatment_comp[meta_mouse_JYNR210602$Treatment_comp == "Chow Vehicle"] <- "Chow Vehicle #1"
meta_mouse_JYNR210602$Treatment_comp[meta_mouse_JYNR210602$Treatment_comp == "DIO-NASH Vehicle"] <- "DIO-NASH Vehicle #1"
meta_mouse_JYNR220601$Treatment_comp[meta_mouse_JYNR220601$Treatment_comp == "Chow Vehicle"] <- "Chow Vehicle #2"
meta_mouse_JYNR220601$Treatment_comp[meta_mouse_JYNR220601$Treatment_comp == "DIO-NASH Vehicle"] <- "DIO-NASH Vehicle #2"

mouse_JYNR210602 <- mouse_JYNR210602[,colnames(mouse_JYNR210602) %in% colnames(mouse_JYNR220601)]
mouse_JYNR220601 <- mouse_JYNR220601[,colnames(mouse_JYNR220601) %in% colnames(mouse_JYNR210602)]
sum(!rownames(mouse_JYNR210602) == rownames(meta_mouse_JYNR210602)) # should be 0
sum(!rownames(mouse_JYNR220601) == rownames(meta_mouse_JYNR220601)) # should be 0
sum(!colnames(mouse_JYNR210602) == colnames(mouse_JYNR220601)) # should be 0

# Create log2-transformed dataframes
log_JYNR210602 <- log2(mouse_JYNR210602)
log_JYNR220601 <- log2(mouse_JYNR220601)

# Create versions without median signal normalization (non-core version of ANML) (and with it)
adat1_norm1 <- mouse_JYNR210602 <- read_adat("~/JDJG/data/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.hybNorm.medNormInt.plateScale.adat")
adat2_norm1 <- mouse_JYNR220601 <- read_adat("~/JDJG/data/raw/Mouse_studies/JYNR220601_GUS2022_898/SS-2231349_v4.1_other.hybNorm.medNormInt.plateScale.adat")
log_1_5_JYNR210602 <- log_JYNR210602
log_1_5_JYNR220601 <- log_JYNR220601

# Align 1-4 dfs with regular dfs
# JYNR210602
cal_JYNR210602 <- mouse_JYNR210602[mouse_JYNR210602$SampleType == "Calibrator",
                                         colnames(mouse_JYNR210602) %in% colnames(log_JYNR210602)]
mouse_JYNR210602 <- mouse_JYNR210602[mouse_JYNR210602$SampleId %in% rownames(log_JYNR210602),]
rownames(mouse_JYNR210602) <- mouse_JYNR210602$SampleId
mouse_JYNR210602 <- mouse_JYNR210602[order(rownames(mouse_JYNR210602)),
                                       colnames(mouse_JYNR210602) %in% colnames(log_JYNR210602)]
log_1_4_JYNR210602 <- log2(mouse_JYNR210602)
sum(!rownames(log_1_4_JYNR210602) == rownames(log_JYNR210602)) # should be 0
sum(!colnames(log_1_4_JYNR210602) == colnames(log_JYNR210602)) # should be 0
old_1_5_JYNR210602 <- log_JYNR210602 # store 1-5 version
log_JYNR210602 <- log_1_4_JYNR210602 # overwrite 1-5 version

# JYNR220601
cal_JYNR220601 <- mouse_JYNR220601[mouse_JYNR220601$SampleType == "Calibrator",
                                         colnames(mouse_JYNR220601) %in% colnames(log_JYNR220601)]
mouse_JYNR220601 <- mouse_JYNR220601[mouse_JYNR220601$SampleId %in% rownames(log_JYNR220601),]
rownames(mouse_JYNR220601) <- mouse_JYNR220601$SampleId
mouse_JYNR220601 <- mouse_JYNR220601[order(rownames(mouse_JYNR220601)),
                                       colnames(mouse_JYNR220601) %in% colnames(log_JYNR220601)]
log_1_4_JYNR220601 <- log2(mouse_JYNR220601)
sum(!rownames(log_1_4_JYNR220601) == rownames(log_JYNR220601)) # should be 0
sum(!colnames(log_1_4_JYNR220601) == colnames(log_JYNR220601)) # should be 0
old_1_5_JYNR220601 <- log_JYNR220601 # store 1-5 version
log_JYNR220601 <- log_1_4_JYNR220601 # overwrite 1-5 version

# Create version with additional outliers removed
out1 <- c(4018696639, 4018711275, 4018700427)
outrm1 <- !meta_mouse_JYNR210602$TubeID %in% out1
outrm_JYNR210602 <- log_JYNR210602[outrm1,]
outrm_meta_JYNR210602 <- meta_mouse_JYNR210602[outrm1,]
out2 <- c(118606)
outrm2 <- !meta_mouse_JYNR220601$`ANIMAL#` %in% out2
outrm_JYNR220601 <- log_JYNR220601[outrm2,]
outrm_meta_JYNR220601 <- meta_mouse_JYNR220601[outrm2,]

## Create version with outliers removed prior to normalisation
# Load normalisation script
source("scripts/soma_norm_function.R")

# Load unnormalised adats
adat1_unnorm <- read_adat("~/JDJG/data/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.adat")
adat2_unnorm <- read_adat("~/JDJG/data/raw/Mouse_studies/JYNR220601_GUS2022_898/SS-2231349_v4.1_other.adat")

# Remove outliers
out1 <- c(4018697918, 4018696639, 4018711275, 4018700427)
keep1 <- !adat1_unnorm$SampleId %in% out1
unnorm1 <- adat1_unnorm[keep1,]; anno1 <- adat1_norm1[keep1,]

out2 <- c(118606)
keep2 <- !adat2_unnorm$SampleId %in% out2
unnorm2 <- adat2_unnorm[keep2,]; anno2 <- adat2_norm1[keep2,]

# Normalisation (steps 1-4)
outrm2_JYNR210602 <- soma_norm(adat = unnorm1,
                               annotation_df = anno1,
                               medNormSMP = FALSE)

outrm2_JYNR220601 <- soma_norm(adat = unnorm2,
                               annotation_df = anno2,
                               medNormSMP = FALSE)

# Wrangle to match current dfs
outrm2_JYNR210602 <- outrm2_JYNR210602[outrm2_JYNR210602$SampleId %in% rownames(outrm_JYNR210602),]
rownames(outrm2_JYNR210602) <- outrm2_JYNR210602$SampleId
outrm2_JYNR210602 <- outrm2_JYNR210602[order(rownames(outrm2_JYNR210602)),
                                       colnames(outrm2_JYNR210602) %in% colnames(outrm_JYNR210602)]
outrm2_JYNR210602 <- log2(outrm2_JYNR210602)

outrm2_JYNR220601 <- outrm2_JYNR220601[outrm2_JYNR220601$SampleId %in% rownames(outrm_JYNR220601),]
rownames(outrm2_JYNR220601) <- outrm2_JYNR220601$SampleId
outrm2_JYNR220601 <- outrm2_JYNR220601[order(rownames(outrm2_JYNR220601)),
                                       colnames(outrm2_JYNR220601) %in% colnames(outrm_JYNR220601)]
outrm2_JYNR220601 <- log2(outrm2_JYNR220601)

# Update df
old_log_JYNR210602 <- log_JYNR210602; old_log_JYNR220601 <- log_JYNR220601
log_JYNR210602 <- outrm2_JYNR210602; log_JYNR220601 <- outrm2_JYNR220601
old_meta_JYNR210602 <- meta_mouse_JYNR210602; old_meta_JYNR220601 <- meta_mouse_JYNR220601
meta_mouse_JYNR210602 <- outrm_meta_JYNR210602; meta_mouse_JYNR220601 <- outrm_meta_JYNR220601

# Get df with no outliers removed
# Normalisation (steps 1-4)
init_JYNR210602 <- adat1_norm1
init_JYNR220601 <- adat2_norm1

# Wrangle to match current dfs
init_JYNR210602 <- init_JYNR210602[init_JYNR210602$SampleType == "Sample",]
rownames(init_JYNR210602) <- init_JYNR210602$SampleId
init_JYNR210602 <- init_JYNR210602[order(rownames(init_JYNR210602)),
                                       colnames(init_JYNR210602) %in% colnames(log_JYNR210602)]
init_JYNR210602 <- log2(init_JYNR210602)

init_JYNR220601 <- init_JYNR220601[init_JYNR220601$SampleType == "Sample",]
rownames(init_JYNR220601) <- init_JYNR220601$SampleId
init_JYNR220601 <- init_JYNR220601[order(rownames(init_JYNR220601)),
                                       colnames(init_JYNR220601) %in% colnames(log_JYNR220601)]
init_JYNR220601 <- log2(init_JYNR220601)

# Create "initial" metadata df
init_meta_JYNR210602 <- old_meta_JYNR210602
init_meta_JYNR210602 <- init_meta_JYNR210602 %>% add_row(TubeID = "4018697918",
                                                         SubjectID = 99099,
                                                         Model = "Healthy",
                                                         Treatment = "Chow Vehicle",
                                                         Group = "1",
                                                         id = "4018697918",
                                                         dataset = "SomaMouse1")
rownames(init_meta_JYNR210602)[48] <- "4018697918"
init_meta_JYNR210602 <- init_meta_JYNR210602[order(rownames(init_meta_JYNR210602)),]
sum(!rownames(init_meta_JYNR210602) == rownames(init_JYNR210602)) # should be 0

init_meta_JYNR220601 <- old_meta_JYNR220601
init_meta_JYNR220601 <- init_meta_JYNR220601 %>% add_row(Treatment = "DIO-NASH Treatment D",
                                                         `ANIMAL#` = "118606",
                                                         Group = "3",
                                                         id = "118606",
                                                         dataset = "SomaMouse2",
                                                         Model = "NASH")
rownames(init_meta_JYNR220601)[41] <- "118606"
init_meta_JYNR220601 <- init_meta_JYNR220601[order(rownames(init_meta_JYNR220601)),]
sum(!rownames(init_meta_JYNR220601) == rownames(init_JYNR220601)) # should be 0

# Get df with no outliers removed
# Normalisation (steps 1-5)
init_1_5_JYNR210602 <- soma_norm(adat = adat1_unnorm,
                             annotation_df = anno1,
                             medNormSMP = TRUE)
init_1_5_JYNR220601 <- soma_norm(adat = adat2_unnorm,
                             annotation_df = anno2,
                             medNormSMP = TRUE)

# Wrangle to match current dfs
init_1_5_JYNR210602 <- init_1_5_JYNR210602[init_1_5_JYNR210602$SampleType == "Sample",]
rownames(init_1_5_JYNR210602) <- init_1_5_JYNR210602$SampleId
init_1_5_JYNR210602 <- init_1_5_JYNR210602[order(rownames(init_1_5_JYNR210602)),
                                   colnames(init_1_5_JYNR210602) %in% colnames(log_JYNR210602)]
init_1_5_JYNR210602 <- log2(init_1_5_JYNR210602)

init_1_5_JYNR220601 <- init_1_5_JYNR220601[init_1_5_JYNR220601$SampleType == "Sample",]
rownames(init_1_5_JYNR220601) <- init_1_5_JYNR220601$SampleId
init_1_5_JYNR220601 <- init_1_5_JYNR220601[order(rownames(init_1_5_JYNR220601)),
                                   colnames(init_1_5_JYNR220601) %in% colnames(log_JYNR220601)]
init_1_5_JYNR220601 <- log2(init_1_5_JYNR220601)

# Get df with steps 1-5
# Normalisation (steps 1-4)
log_1_5_JYNR210602 <- soma_norm(adat = unnorm1,
                               annotation_df = anno1,
                               medNormSMP = TRUE)

log_1_5_JYNR220601 <- soma_norm(adat = unnorm2,
                               annotation_df = anno2,
                               medNormSMP = TRUE)

# Wrangle to match current dfs
log_1_5_JYNR210602 <- log_1_5_JYNR210602[log_1_5_JYNR210602$SampleId %in% rownames(outrm_JYNR210602),]
rownames(log_1_5_JYNR210602) <- log_1_5_JYNR210602$SampleId
log_1_5_JYNR210602 <- log_1_5_JYNR210602[order(rownames(log_1_5_JYNR210602)),
                                       colnames(log_1_5_JYNR210602) %in% colnames(outrm_JYNR210602)]
log_1_5_JYNR210602 <- log2(log_1_5_JYNR210602)

log_1_5_JYNR220601 <- log_1_5_JYNR220601[log_1_5_JYNR220601$SampleId %in% rownames(outrm_JYNR220601),]
rownames(log_1_5_JYNR220601) <- log_1_5_JYNR220601$SampleId
log_1_5_JYNR220601 <- log_1_5_JYNR220601[order(rownames(log_1_5_JYNR220601)),
                                       colnames(log_1_5_JYNR220601) %in% colnames(outrm_JYNR220601)]
log_1_5_JYNR220601 <- log2(log_1_5_JYNR220601)

# Perform cross-study normalization
source("~/JDJG/scripts/wrangle_crossnorm.R")

# # Scale all aptamers to be within (0-10)
# # (and make sure the scale factors are the same across datasets)
# source("~/JDJG/scripts/scale_function.R")
# # scale all the dfs
# scaled_log_JYNR210602 <- scale_function(log_JYNR210602, log_JYNR220601)
# scaled_log_JYNR220601 <- scale_function(log_JYNR220601, log_JYNR210602)
# scaled_norm_JYNR210602 <- scale_function(norm_JYNR210602, norm_JYNR220601)
# scaled_norm_JYNR220601 <- scale_function(norm_JYNR220601, norm_JYNR210602)
# scaled_limma_JYNR210602 <- scale_function(limma_JYNR210602, limma_JYNR220601)
# scaled_limma_JYNR220601 <- scale_function(limma_JYNR220601, limma_JYNR210602)

# Compute Z-scores
# Super slow to run
# For each SOMAmer:
mouse_JYNR210602 <- mouse_JYNR210602[rownames(mouse_JYNR210602) %in% rownames(log_JYNR210602),]
mouse_JYNR220601 <- mouse_JYNR220601[rownames(mouse_JYNR220601) %in% rownames(log_JYNR220601),]
# z_JYNR210602 <- mouse_JYNR210602; z_JYNR220601 <- mouse_JYNR220601
# for (i in 1:ncol(mouse_JYNR210602)){
#   print(i)
#   mean_JYNR210602 <- mean(mouse_JYNR210602[,i])
#   sd_JYNR210602 <- sd(mouse_JYNR210602[,i])
#   mean_JYNR220601 <- mean(mouse_JYNR220601[,i])
#   sd_JYNR220601 <- sd(mouse_JYNR220601[,i])
#   for (j in 1:nrow(mouse_JYNR210602)){
#     z_JYNR210602[j,i] <- (mouse_JYNR210602[j,i] - mean_JYNR210602) / sd_JYNR210602
#   }
#   for (j in 1:nrow(mouse_JYNR220601)){
#     z_JYNR220601[j,i] <- (mouse_JYNR220601[j,i] - mean_JYNR220601) / sd_JYNR220601
#   }
# }
# 
# # Compute Z-scores (based on log2-transformed data)
# # Super slow to run
# # For each SOMAmer:
# z_log_JYNR210602 <- log_JYNR210602; z_log_JYNR220601 <- log_JYNR220601
# for (i in 1:ncol(log_JYNR210602)){
#   print(i)
#   mean_JYNR210602 <- mean(log_JYNR210602[,i])
#   sd_JYNR210602 <- sd(log_JYNR210602[,i])
#   mean_JYNR220601 <- mean(log_JYNR220601[,i])
#   sd_JYNR220601 <- sd(log_JYNR220601[,i])
#   for (j in 1:nrow(log_JYNR210602)){
#     z_log_JYNR210602[j,i] <- (log_JYNR210602[j,i] - mean_JYNR210602) / sd_JYNR210602
#   }
#   for (j in 1:nrow(log_JYNR220601)){
#     z_log_JYNR220601[j,i] <- (log_JYNR220601[j,i] - mean_JYNR220601) / sd_JYNR220601
#   }
# }
# 
# # # Z-scores across datasets (matches what's being done in the PCA)
# z_joint_JYNR210602 <- log_JYNR210602; z_joint_JYNR220601 <- log_JYNR220601
# for (i in 1:ncol(log_JYNR210602)){
#   print(i)
#   mean_joint <- mean(c(log_JYNR210602[,i], log_JYNR220601[,i]))
#   sd__joint <- sd(c(log_JYNR210602[,i], log_JYNR220601[,i]))
#   for (j in 1:nrow(log_JYNR210602)){
#     z_joint_JYNR210602[j,i] <- (log_JYNR210602[j,i] - mean_joint) / sd__joint
#   }
#   for (j in 1:nrow(log_JYNR220601)){
#     z_joint_JYNR220601[j,i] <- (log_JYNR220601[j,i] - mean_joint) / sd__joint
#   }
# }

## Save
save("meta_mouse_JYNR210602", "old_meta_JYNR210602", "init_meta_JYNR210602",
     "init_JYNR210602", "old_log_JYNR210602", "log_JYNR210602",
     "init_1_5_JYNR210602", "old_1_5_JYNR210602", "log_1_5_JYNR210602",
     "crossnorm_JYNR210602", "limma_JYNR210602",
     "z_joint_JYNR210602",
     "cal_one_JYNR210602", "cal_three_JYNR210602", "cal_indi_JYNR210602",
     "norm_JYNR210602", "z_JYNR210602", "z_log_JYNR210602",
     
     "meta_mouse_JYNR220601","old_meta_JYNR220601", "init_meta_JYNR220601",
     "init_JYNR220601", "old_log_JYNR220601", "log_JYNR220601",
     "init_1_5_JYNR220601", "old_1_5_JYNR220601", "log_1_5_JYNR220601",
     "crossnorm_JYNR220601", "limma_JYNR220601",
     "z_joint_JYNR220601",
     "norm_JYNR220601", "z_JYNR220601", "z_log_JYNR220601",
     
     "SomaScanAnnotation", "case_control_norm_dfs",
     file = "~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")

# Clean up
rm(list = c("Flag_n_dropListInfo_Mouse_JYNR210602",
            "Flag_n_dropListInfo_Mouse_JYNR220601",
            "Phenotype_Mouse_JYNR210602",
            "Phenotype_Mouse_JYNR210602_w_soma",
            "Phenotype_Mouse_JYNR220601",
            "Phenotype_Mouse_JYNR220601_w_soma",
            "mouse_JYNR210602", "mouse_JYNR220601",
            "log_1_4_JYNR210602", "log_1_4_JYNR220601",
            "cal_JYNR210602", "cal_JYNR220601",
            "adat1_unnorm", "adat1_norm1", "anno1", "out1", "keep1",
            "adat2_unnorm", "adat2_norm1", "anno2", "out2", "keep2",
            "outrm_meta_JYNR210602", "outrm_meta_JYNR220601",
            "outrm2_JYNR210602", "outrm2_JYNR220601"))


# ################################################################################
# ### Minipig data
# 
# ## Load data
# # load("~/NNEDL/masterprojectxjdjg/curated/MiniPig_studies/curated_data_minipigs.RData")
# load("~/JDJG/data/curated/MiniPig_studies/curated_data_minipigs.RData")
# 
# # Get annotation file
# SomaScanAnnotation <- read.csv("../NNEDL/somalogicannotation/curated/v4_1_7k/current/SomaAnnotation_v4_1_7k_2023-11-14.csv")
# SomaScanAnnotation$SeqId <- paste0("seq.",sub("\\-.*", "", SomaScanAnnotation$SeqId),".",sub('.+-(.+)', '\\1', SomaScanAnnotation$SeqId))
# 
# ## Wrangle
# # BQC210206 (time series)
# meta_pig_BQC210206 <- Phenotype_Minipig_BQC230103 %>% 
#   filter(Unique_ID %in% Phenotype_Minipig_BQC230103_w_soma$Unique_ID)
# meta_pig_BQC210206$id <- meta_pig_BQC210206$Unique_ID
# meta_pig_BQC210206 <- meta_pig_BQC210206[order(meta_pig_BQC210206$Unique_ID),]
# meta_pig_BQC210206 <- meta_pig_BQC210206[order(as.numeric(meta_pig_BQC210206$Time_day)),]
# meta_pig_BQC210206 <- meta_pig_BQC210206[order(meta_pig_BQC210206$Fasting),]
# rownames(meta_pig_BQC210206) <- meta_pig_BQC210206$Unique_ID
# pig_BQC210206 <- Phenotype_Minipig_BQC230103_w_soma
# pig_BQC210206 <- pig_BQC210206[order(pig_BQC210206$Unique_ID),]
# pig_BQC210206 <- pig_BQC210206[order(as.numeric(pig_BQC210206$Time_day)),]
# pig_BQC210206 <- pig_BQC210206[order(pig_BQC210206$Fasting),]
# rownames(pig_BQC210206) <- pig_BQC210206$Unique_ID
# pig_BQC210206 <- pig_BQC210206[,-c(1, 7258:7270)]
# dat <- pig_BQC210206
# for (i in 1:dim(dat)[2]){
#   if (!(str_detect(colnames(dat)[i], "seq"))){
#     print(i)
#     print(colnames(dat)[i])
#   }
# }
# # Make "Group" column with number IDs
# colnames(meta_pig_BQC210206)[9] <- "Treatment" # rename "Group" to "Treatment"
# meta_pig_BQC210206 <- meta_pig_BQC210206 %>% 
#   mutate(Group = case_when(
#     Time_day_fac == -3 & Treatment == "Vehicle" ~ "1a",
#     Time_day_fac == -3 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "1b",
#     Time_day_fac == 22 & Treatment == "Vehicle" ~ "2a",
#     Time_day_fac == 22 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "2b",
#     Time_day_fac == 50 & Treatment == "Vehicle" ~ "3a",
#     Time_day_fac == 50 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "3b",
#     Time_day_fac == 73 & Treatment == "Vehicle" ~ "4a",
#     Time_day_fac == 73 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "4b",
#     Time_day_fac == 86 & Treatment == "Vehicle" ~ "5a",
#     Time_day_fac == 86 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "5b",
#     Time_day_fac == 1 & Treatment == "Vehicle" ~ "6a",
#     Time_day_fac == 1 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "6b",
#     Time_day_fac == 3 & Treatment == "Vehicle" ~ "7a",
#     Time_day_fac == 3 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "7b",
#     Time_day_fac == 10 & Treatment == "Vehicle" ~ "8a",
#     Time_day_fac == 10 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "8b",
#     Time_day_fac == 24 & Treatment == "Vehicle" ~ "9a",
#     Time_day_fac == 24 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "9b",
#     Time_day_fac == 38 & Treatment == "Vehicle" ~ "10a",
#     Time_day_fac == 38 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "10b",
#     Time_day_fac == 59 & Treatment == "Vehicle" ~ "11a",
#     Time_day_fac == 59 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "11b",
#     Time_day_fac == 81 & Treatment == "Vehicle" ~ "12a",
#     Time_day_fac == 81 & Treatment == "32.5 nmol/kg - 0518-9336-0768-5K" ~ "12b"
#   ))
# meta_pig_BQC210206$Group <- as_factor(meta_pig_BQC210206$Group)
# 
# # Fix animal IDs
# meta_pig_BQC210206$Animal_ID[meta_pig_BQC210206$Animal_ID == "935"] <- "235-935"
# meta_pig_BQC210206$Animal_ID[meta_pig_BQC210206$Animal_ID == "342-987"] <- "235-987"
# 
# 
# # BQC230103 (tox)
# meta_pig_BQC230103 <- as.data.frame(Phenotype_Minipig_Tox) %>% 
#   filter(Unique_ID %in% Phenotype_Minipig_Tox_w_soma$Unique_ID)
# meta_pig_BQC230103$id <- meta_pig_BQC230103$Unique_ID
# meta_pig_BQC230103 <- meta_pig_BQC230103[order(meta_pig_BQC230103$Short_ID),]
# meta_pig_BQC230103 <- meta_pig_BQC230103[order(meta_pig_BQC230103$`Time Point`),]
# rownames(meta_pig_BQC230103) <- meta_pig_BQC230103$id
# pig_BQC230103 <- Phenotype_Minipig_Tox_w_soma
# pig_BQC230103 <- pig_BQC230103[order(pig_BQC230103$Short_ID),]
# pig_BQC230103 <- pig_BQC230103[order(pig_BQC230103$`Time Point`),]
# sum(!rownames(meta_pig_BQC230103) == rownames(pig_BQC230103)) # should be 0
# pig_BQC230103 <- pig_BQC230103[,-c(1,7258:7285)]
# dat <- pig_BQC230103
# for (i in 1:dim(dat)[2]){
#   if (!(str_detect(colnames(dat)[i], "seq"))){
#     print(i)
#     print(colnames(dat)[i])
#   }
# }
# 
# # Make "Group" column with number IDs
# colnames(meta_pig_BQC230103)[14] <- "Treatment" # rename "Group" to "Treatment"
# meta_pig_BQC230103 <- meta_pig_BQC230103 %>% 
#   mutate(Treatment = case_when(
#     Treatment == 1 ~ "Vehicle",
#     Treatment == 3 ~ "2 mg/kg",
#     Treatment == 4 ~ "8 mg/kg"
#     ))
# meta_pig_BQC230103 <- meta_pig_BQC230103 %>% 
#   mutate(Group = case_when(
#     Time_day_fac == 3 & Treatment == "Vehicle" ~ "1a",
#     Time_day_fac == 3 & Treatment == "2 mg/kg" ~ "1b",
#     Time_day_fac == 3 & Treatment == "8 mg/kg" ~ "1c",
#     Time_day_fac == 9 & Treatment == "Vehicle" ~ "2a",
#     Time_day_fac == 9 & Treatment == "2 mg/kg" ~ "2b",
#     Time_day_fac == 9 & Treatment == "8 mg/kg" ~ "2c",
#     Time_day_fac == 93 & Treatment == "Vehicle" ~ "3a",
#     Time_day_fac == 93 & Treatment == "2 mg/kg" ~ "3b",
#     Time_day_fac == 93 & Treatment == "8 mg/kg" ~ "3c",
#   ))
# meta_pig_BQC230103$Group <- as_factor(meta_pig_BQC230103$Group)
# 
# 
# # Align dataframes
# pig_BQC210206 <- pig_BQC210206[,colnames(pig_BQC210206) %in% colnames(pig_BQC230103)]
# pig_BQC230103 <- pig_BQC230103[,colnames(pig_BQC230103) %in% colnames(pig_BQC210206)]
# sum(!rownames(pig_BQC210206) == rownames(meta_pig_BQC210206)) # should be 0
# rownames(meta_pig_BQC230103) <- meta_pig_BQC230103$id
# sum(!rownames(pig_BQC230103) == rownames(meta_pig_BQC230103)) # should be 0
# sum(!colnames(pig_BQC210206) == colnames(pig_BQC230103)) # should be 0
# 
# # Create log2-transformed dataframes
# log_BQC210206 <- log2(pig_BQC210206)
# log_BQC230103 <- log2(pig_BQC230103)
# 
# ## Save
# save("meta_pig_BQC210206",
#      "pig_BQC210206", "log_BQC210206",
#      "meta_pig_BQC230103",
#      "pig_BQC230103", "log_BQC230103",
#      "SomaScanAnnotation",
#      file = "~/JDJG/data/curated/MiniPig_studies/curated_data_minipigs_2.RData")
