## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel)
library(pheatmap); library(limma)
setwd("~/JDJG/clean/")
source("scripts/SomaScan_Mouse/overlap_function.R")

## Load data and perform initial wrangling (this should be further automated for future use)
# load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")


################################################################################
### Analysis
# soma_data_1 <- z_JYNR210602
# soma_data_2 <- z_JYNR220601
# heatmap_data_1 <- t(as.matrix(soma_data_1[meta_mouse_JYNR210602$Group %in% c(3,1),]))
# 
# 
# heatmap_data_1 <- t(as.matrix(soma_data_1))
# heatmap_data_2 <- t(as.matrix(soma_data_2[meta_mouse_JYNR220601$Group %in% c(2,1),]))
# 
# heatmap(heatmap_data_1, ColSideColors = meta_mouse_JYNR210602$Group)
# heatmap(heatmap_data_2)

# # Heatmap
# # Individual datasets
# for (dataset in (c("JYNR210602", "JYNR220601"))){
#   # do stuff
#   for (handle in c("log_", "z_", "z_log_", "z_joint_")){
#     # do stuff
#     soma_data <- t(as.matrix(get(paste0(handle, dataset))))
#     sampleinfo <- get(paste0("meta_mouse_", dataset)) %>%
#       select(Model, Group)
#     corMatrix <- cor(soma_data,use = "c")
#     rownames(sampleinfo) <- colnames(corMatrix)
#     # identical(rownames(sampleinfo1), colnames(corMatrix1))
#     if (dataset == "JYNR210602"){
#       groups <- c("1" = "#74DFE2", "2" = "#20979B",
#                   "3" = "#F7A8A3", "4" = "#D0756F", "5" = "#A9433C", "6" = "#821008")
#     } else {
#       groups <- c("1" = "#74DFE2", "2" = "#F7A8A3", "3" = "#D2B8B5"
#                   )
#     }
#     ann_colours <- list(
#       Model = c(Healthy = "#2EBABF", NASH = "#BF3D35"),
#       Group = groups)
#     if (handle == "log_"){
#       title <- paste0(dataset, ", log2-transformed RFU values")
#     }
#     if (handle == "z_"){
#       title <- paste0(dataset, ", z-scores")
#     }
#     if (handle == "z_log_"){
#       title <- paste0(dataset, ", z-scores (based on log2-transformed data)")
#     }
#     if (handle == "z_joint_"){
#       title <- paste0(dataset, ", z-scores (across datasets, log2 data)")
#     }
#     pheatmap(corMatrix, annotation_col = sampleinfo,
#              show_rownames = FALSE, show_colnames = FALSE, annotation_colors = ann_colours,
#              filename = paste0("figures/heatmap_", handle, dataset, ".png"), main = title,
#              width = 10.5, height = 11.5, fontsize = 16)
#   }
# }

# Heatmap
# Both datasets together
# # for (handle in c("log_", "z_", "z_log_", "z_joint_")){
# all_groups <- FALSE
# all_groups <- TRUE
# for (handle in c("init_", "log_", "init_1_5_", "log_1_5_", "z_joint_")){
# # for (handle in c("z_joint_")){
#   # for (correction in c("", "_limma")){
#   for (correction in c("")){
#     soma_data_1 <- t(as.matrix(get(paste0(handle, "JYNR210602"))))
#     soma_data_2 <- t(as.matrix(get(paste0(handle, "JYNR220601"))))
#     soma_data_merged <- merge(soma_data_1, soma_data_2, by = "row.names", all = TRUE)
#     rownames(soma_data_merged) <- soma_data_merged[,1]
#     soma_data_merged <- soma_data_merged[-1]
#     batch <- c(rep("Data #1",dim(soma_data_1)[2]),rep("Data #2",dim(soma_data_2)[2]))
#     if (correction == "_limma"){
#       soma_data_merged <- removeBatchEffect(soma_data_merged,batch)
#     }
#     
#     # Get metadata
#     if (handle %in% c("init_", "init_1_5_")){
#       meta1 <- init_meta_JYNR210602; meta2 <- init_meta_JYNR220601
#     } else if (handle %in% c("log_", "log_1_5_", "z_joint_")){
#       meta1 <- meta_mouse_JYNR210602; meta2 <- meta_mouse_JYNR220601
#     } else {
#       print("METADATA NOT IDENTIFIED BY HANDLE !")
#       break
#     }
#     
#     if (all_groups == TRUE){
#       meta1_cc <- meta1[,colnames(meta1) %in% colnames(meta2)]
#       meta2_cc <- meta2[,colnames(meta2) %in% colnames(meta1)]
#       grps <- c(rownames(meta1), rownames(meta2))
#     } else {
#       meta1_cc <- meta1[meta1$Group %in% c(1,3), colnames(meta1) %in% colnames(meta2)]
#       meta2_cc <- meta2[meta2$Group %in% c(1,2), colnames(meta2) %in% colnames(meta1)]
#       grps <- c(rownames(meta1[meta1$Group %in% c(1,3),]),
#                 rownames(meta2[meta2$Group %in% c(1,2),]))
#     }
#     # Keep only case/control vehicle groups
#     soma_data_merged <- soma_data_merged[,colnames(soma_data_merged) %in% grps]
#     meta_merged <- dplyr::full_join(meta1_cc, meta2_cc)
#     colnames(meta_merged)[colnames(meta_merged) == "dataset"] <- "Dataset"
#     # meta_merged$Group[meta_merged$id == "118499"] <- 9 # comment out
#     
#     sampleinfo <- meta_merged %>%
#       select(Model, Dataset)
#     corMatrix <- cor(soma_data_merged,use = "c")
#     rownames(sampleinfo) <- colnames(corMatrix) 
#     
#     
#     ann_colours <- list(Dataset = c(SomaMouse1 = "#FFA355",
#                                     SomaMouse2 = "#55B1FF"),
#                         Model = c(Healthy = "#53BF2E", NASH = "#BF3D35"))
#     if (handle == "log_"){
#       title <- paste0("Outliers removed")
#     }
#     if (handle == "init_"){
#       title <- paste0("No outliers removed")
#     }
#     if (handle == "log_1_5_"){
#       title <- paste0("Outliers removed, median normalised")
#     }
#     if (handle == "init_1_5_"){
#       title <- paste0("No outliers removed, median normalised")
#     }
#     
#     if (handle == "z_"){
#       title <- paste0("Z-scores")
#     }
#     if (handle == "z_log_"){
#       title <- paste0("Z-scores (log2)")
#     }
#     if (handle == "z_joint_"){
#       title <- paste0("Outliers removed, Z-scores")
#     }
#     if (correction == "_limma"){
#       title <- paste0(title, ", limma")
#     }
#     filename <- paste0("figures/heatmap_", handle, "both_datasets", correction)
#     if (all_groups == TRUE){
#       filename <- paste0(filename, "_allgrps")
#     }
#     filename <- paste0(filename, ".png")
#     pheatmap(corMatrix, annotation_col = sampleinfo,
#              show_rownames = FALSE, show_colnames = FALSE, annotation_colors = ann_colours,
#              filename = filename,
#              main = title, width = 10.5, height = 11.5, fontsize = 16)
#   }
# }



# # Heatmap
# # Both datasets together, all groups
# for (handle in c("log_", "z_joint_")){
#   # for (handle in c("z_joint_")){
#   for (correction in c("", "_limma")){
#     soma_data_1 <- t(as.matrix(get(paste0(handle, "JYNR210602"))))
#     soma_data_2 <- t(as.matrix(get(paste0(handle, "JYNR220601"))))
#     soma_data_merged <- merge(soma_data_1, soma_data_2, by = "row.names", all = TRUE)
#     rownames(soma_data_merged) <- soma_data_merged[,1]
#     soma_data_merged <- soma_data_merged[-1]
#     batch <- c(rep("Data #1",dim(soma_data_1)[2]),rep("Data #2",dim(soma_data_2)[2]))
#     if (correction == "_limma"){
#       soma_data_merged <- removeBatchEffect(soma_data_merged,batch)
#     }
#     
#     # Keep only case/control vehicle groups
#     grps <- c(rownames(meta_mouse_JYNR210602),
#               rownames(meta_mouse_JYNR220601))
#     soma_data_merged <- soma_data_merged[,colnames(soma_data_merged) %in% grps]
#     
#     meta1 <- meta_mouse_JYNR210602[,
#                                    colnames(meta_mouse_JYNR210602) %in% colnames(meta_mouse_JYNR220601)]
#     meta2 <- meta_mouse_JYNR220601[,
#                                    colnames(meta_mouse_JYNR220601) %in% colnames(meta_mouse_JYNR210602)]
#     meta_merged <- dplyr::full_join(meta1, meta2)
#     meta_merged$outlier <- 0
#     meta_merged$outlier[meta_merged$id == "118499"] <- 1 # comment out
#     
#     sampleinfo <- meta_merged %>%
#       select(Model, dataset, Group, outlier) # change these, as "Group" overlaps (but the 2s and 3s are not equivalent across studies)
#     corMatrix <- cor(soma_data_merged,use = "c")
#     rownames(sampleinfo) <- colnames(corMatrix) 
#     
#     
#     ann_colours <- list(
#       Model = c(Healthy = "#2EBABF", NASH = "#BF3D35"))
#     if (handle == "log_"){
#       title <- paste0("Log2-transformed")
#     }
#     if (handle == "z_"){
#       title <- paste0("Z-scores")
#     }
#     if (handle == "z_log_"){
#       title <- paste0("Z-scores (log2)")
#     }
#     if (handle == "z_joint_"){
#       title <- paste0("Z-scores (across, log2)")
#     }
#     if (correction == "_limma"){
#       title <- paste0(title, ", limma")
#     }
#     pheatmap(corMatrix, annotation_col = sampleinfo,
#              show_rownames = FALSE, show_colnames = FALSE, annotation_colors = ann_colours,
#              filename = paste0("figures/heatmap_", handle, "both_datasets", correction, "_all_grps", ".png"),
#              main = title, width = 10.5, height = 11.5, fontsize = 16)
#   }
# }
