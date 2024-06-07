### Initialize document
# Set working directory
setwd("~/JDJG/clean/")


# Load functions
source("scripts/SomaScan_Mouse/pca_function.R")


### LOAD DATA
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")


### Pre-process data
meta1 <- init_meta_JYNR210602[,colnames(init_meta_JYNR210602) %in% colnames(init_meta_JYNR220601)]
meta2 <- init_meta_JYNR220601[,colnames(init_meta_JYNR220601) %in% colnames(init_meta_JYNR210602)]
comb_meta <- dplyr::full_join(meta1, meta2)
comb_meta$label <- paste(comb_meta$dataset, comb_meta$Treatment_comp, sep = ", ")
colnames(comb_meta)[colnames(comb_meta) == "dataset"] <- "Dataset" # rename "dataset" to "Dataset"
comb_meta$Treatment <- factor(comb_meta$Treatment,
                              levels = c("Chow Vehicle", "Chow Treatment A",
                                         "DIO-NASH Vehicle", "DIO-NASH Treatment A",
                                         "DIO-NASH Treatment B", "DIO-NASH Treatment C",
                                         "DIO-NASH Treatment D"))
old_outrm <- comb_meta$id %in% c(old_meta_JYNR210602$id, old_meta_JYNR220601$id)
new_outrm <- comb_meta$id %in% c(meta_mouse_JYNR210602$id, meta_mouse_JYNR220601$id)
comb_meta_init <- comb_meta
comb_meta_old <- comb_meta[old_outrm,]
comb_meta <- comb_meta[new_outrm,]
comb_init <- merge_dfs(init_JYNR210602, init_JYNR220601)
comb_old <- merge_dfs(old_log_JYNR210602, old_log_JYNR220601)
comb_baseline <- merge_dfs(log_JYNR210602, log_JYNR220601)
comb_init_1_5 <- merge_dfs(init_1_5_JYNR210602, init_1_5_JYNR220601)
comb_old_1_5 <- merge_dfs(old_1_5_JYNR210602, old_1_5_JYNR220601)
comb_1_5 <- merge_dfs(log_1_5_JYNR210602, log_1_5_JYNR220601)
comb_crossnorm <- merge_dfs(crossnorm_JYNR210602, crossnorm_JYNR220601)
comb_limma <- merge_dfs(limma_JYNR210602, limma_JYNR220601)

# Extra (z-scores and z-scores (log2))
comb_z <- merge_dfs(z_JYNR210602, z_JYNR220601)
comb_z_log <- merge_dfs(z_log_JYNR210602, z_log_JYNR220601)
comb_z_joint <- merge_dfs(z_joint_JYNR210602, z_joint_JYNR220601)

## Create lists with both RFU and meta data
init_data <- list(comb_init, comb_meta_init)
old_data <- list(comb_old, comb_meta_old)
baseline_data <- list(comb_baseline, comb_meta)
init_1_5_data <- list(comb_init_1_5, comb_meta_init)
old_1_5_data <- list(comb_old_1_5, comb_meta_old)
log_1_5_data <- list(comb_1_5, comb_meta)
crossnorm_data <- list(comb_crossnorm, comb_meta)
limma_data <- list(comb_limma, comb_meta)
z_data <- list(comb_z, comb_meta); z_log_data <- list(comb_z_log, comb_meta)
z_joint_data <- list(comb_z_joint, comb_meta)

### PCA
# ## Initial PCAs
# # Baseline and outlier removal
# pcaplot(init_data, title = "Mouse data, no outliers removed", legend = "left", colour = "Treatment")
# ggsave("figures/pca_init.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# pcaplot(old_data, title = "Mouse data, two outliers removed", legend = "left", colour = "Treatment")
# ggsave("figures/pca_old.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# pcaplot(baseline_data, title = "Mouse data, five outliers removed", legend = "left", colour = "Treatment")
# ggsave("figures/pca_base_initial.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# 
# # Baseline and outlier removal (steps 1-5)
# pcaplot(init_1_5_data, title = "Mouse data, no outliers removed, interplate median signal normalisation", legend = "left", colour = "Treatment")
# ggsave("figures/pca_1_5_init.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# pcaplot(old_1_5_data, title = "Mouse data, two outliers removed, interplate median signal normalisation", legend = "left", colour = "Treatment")
# ggsave("figures/pca_1_5_old.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# pcaplot(log_1_5_data, title = "Mouse data, five outliers removed, interplate median signal normalisation", legend = "left", colour = "Treatment")
# ggsave("figures/pca_1_5.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")


## Thesis version of the plot (fig. 1B)
myplots <- list(pcaplot(init_data, title = "No outliers removed",
                        legend = "none", colour = "Treatment", fig1 = TRUE),
                pcaplot(init_1_5_data, title = "No outliers removed, interplate median signal normalisation",
                        legend = "none", colour = "Treatment", fig1 = TRUE),
                pcaplot(baseline_data, title = "Outliers removed",
                        legend = "none", colour = "Treatment", fig1 = TRUE))
legend <- pcaplot(init_data, colour = "Treatment", legend = "top", get_legend = TRUE)
pca_final <- ggarrange(plotlist = myplots,
                       nrow = length(myplots),
                       common.legend = TRUE,
                       legend.grob = legend,
                       legend = "left",
                       hjust = -2)
# pca_final
# ggsave("figures/pca_final.png", width = 3000, height = 2000, units = "px", dpi = 200, bg = "white")

pca_cv_final <- ggarrange(plotlist = list(pca_final, cv_final),
          widths = c(3,1))
# pca_cv_final
# ggsave("figures/pca_cv_final.png", width = 4000, height = 2000, units = "px", dpi = 200, bg = "white")
# 4000 x 2000 -> 3000 x 2000


# ## Thesis version of the plot (fig. 3C)
# myplots <- list(pcaplot(baseline_data, title = "Baseline",
#                         legend = "none", colour = "Treatment"),
#                 pcaplot(crossnorm_data, title = "Cross-normalisation (SomaLogic)",
#                         legend = "none", colour = "Treatment"),
#                 pcaplot(limma_data, title = "Batch effect correction (limma)",
#                         legend = "none", colour = "Treatment"))
# legend <- pcaplot(baseline_data, colour = "Treatment", legend = "top", get_legend = TRUE)
# pca_sup <- ggarrange(plotlist = myplots,
#                        nrow = length(myplots),
#                        common.legend = TRUE,
#                        legend.grob = legend,
#                        legend = "left",
#                        hjust = -2)
# pca_sup


# # All groups
# # pcaplot(log_data, title = "Mouse data, baseline", legend = "left", colour = "Treatment")
# # ggsave("figures/pca_log_treat_hires.png", width = 3000, height = 1600, units = "px", dpi = 200, bg = "white")
# # ggsave("figures/pca_log_treat.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# # pcaplot(crossnorm_data, title = "Mouse data, cross-normalized (SomaLogic)", legend = "left", colour = "Treatment")
# # ggsave("figures/pca_crossnorm_treat.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# # pcaplot(limma_data, title = "Mouse data, batch effects corrected (limma)", legend = "left", colour = "Treatment")
# # ggsave("figures/pca_limma_treat.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# 
# # Only overlapping groups, study/model
# pcaplot(log_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), title = "Mouse data, baseline, only overlapping control groups", legend = "left")
# ggsave("figures/pca_log_overlap_hires.png", width = 3000, height = 1600, units = "px", dpi = 200, bg = "white")
# ggsave("figures/pca_log_overlap.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# pcaplot(crossnorm_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), title = "Mouse data, cross-normalized (SomaLogic), only overlapping control groups", legend = "left")
# ggsave("figures/pca_crossnorm_overlap.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# pcaplot(limma_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), title = "Mouse data, batch effects corrected (limma), only overlapping control groups", legend = "left")
# ggsave("figures/pca_limma_overlap.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# 
# 
# # 1-5 data
# pcaplot(log_1_5_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), title = "Mouse data, baseline, steps 1-5, only overlapping control groups", legend = "left")
# ggsave("figures/pca_log_1_5_overlap.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# pcaplot(log_1_5_data, title = "Mouse data, baseline, steps 1-5", legend = "left", colour = "Treatment")
# ggsave("figures/pca_log_1_5_treat.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# 
# # Outlier removal data
# pcaplot(outrm_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), title = "Mouse data, baseline, outliers removed, only overlapping control groups", legend = "left")
# ggsave("figures/pca_outrm_overlap.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# pcaplot(outrm_data, title = "Mouse data, baseline, outliers removed", legend = "left", colour = "Treatment")
# ggsave("figures/pca_outrm_treat.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# pcaplot(outrm2_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), title = "Mouse data, baseline, outliers removed pre-norm, only overlapping control groups", legend = "left")
# ggsave("figures/pca_outrm2_overlap.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# pcaplot(outrm2_data, title = "Mouse data, baseline, outliers removed pre-norm", legend = "left", colour = "Treatment")
# ggsave("figures/pca_outrm2_treat.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# 
# 
# # # Z-scores
# # pcaplot(z_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), title = "Mouse data, z-scores, only overlapping control groups", legend = "left")
# # ggsave("figures/pca_z_overlap.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# # pcaplot(z_log_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), title = "Mouse data, z-scores (log2), only overlapping control groups", legend = "left")
# # ggsave("figures/pca_z_log_overlap.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# # pcaplot(z_joint_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), title = "Mouse data, z-scores across datasets (log2), only overlapping control groups", legend = "left")
# # ggsave("figures/pca_z_joint_overlap.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# # pcaplot(z_joint_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), title = "Mouse data, z-scores across datasets (log2), only overlapping control groups, no scaling or centering", legend = "left",
# #         scale = FALSE, center = FALSE)
# # ggsave("figures/pca_z_joint_overlap_no_scaling.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
#   # 
#   # # All treatment groups (to show effect/non-effect of z-scores/PCA)
#   # pcaplot(z_joint_data, title = "Mouse data, z-scores across datasets (log2)", colour = "Treatment", legend = "left")
#   # ggsave("figures/pca_z_joint_treat.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
#   # pcaplot(z_joint_data, title = "Mouse data, z-scores across datasets (log2), no scaling or centering",
#   #         colour = "Treatment", legend = "left",
#   #         scale = FALSE, center = FALSE)
#   # ggsave("figures/pca_z_joint_treat_no_scaling.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# # 
# Individual datasets, treatment
# pcaplot(list(log_JYNR210602, comb_meta[comb_meta$Dataset == "SomaMouse1",]),
#         title = "SomaMouse1, baseline",
#         legend = "left", colour = "Treatment",
#         scale = TRUE, center = TRUE)
# ggsave("figures/pca_JYNR210602_treat.png", width = 3000, height = 800, units = "px", dpi = 200, bg = "white")
# pcaplot(list(log_JYNR220601, comb_meta[comb_meta$Dataset == "SomaMouse2",]),
#         title = "SomaMouse2, baseline",
#         legend = "left", colour = "Treatment",
#         scale = TRUE, center = TRUE)
# ggsave("figures/pca_JYNR220601_treat.png", width = 3000, height = 800, units = "px", dpi = 200, bg = "white")

# Fuck around and find out
# pcaplot(z_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"),
#         title = "Mouse data, z-scores, only overlapping control groups",
#         colour = "Treatment_comp", legend = "left",
#         scale = TRUE, center = TRUE)
# pcaplot(log_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"),
#         title = "Mouse data, baseline, only overlapping control groups",
#         legend = "left",
#         scale = TRUE, center = TRUE)
#
# pcaplot(log_data, title = "Mouse data, baseline", colour = "Treatment_comp", legend = "left")
# ggsave("figures/pca_log_treat.png", width = 3000, height = 1600, units = "px", dpi = 200, bg = "white")
#
#
# pcaplot(log_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"), colour = "dataset", title = "Mouse data, baseline, only overlapping control groups", legend = "left")
# pcaplot(log_data)
# 
# 
# 
# 
# 
# 
# 
# 
# ################################################################################
# ## PCA - protein weights
# # For baseline with outliers removed (pre-norm), PC4 accurately separates the two datasets
# # Hence, we want to look at the proteins with the largest weights on that axis
# # These are the ones driving the inter-study differences.
# soma <- outrm2_data[[1]]; meta <- outrm2_data[[2]]
# soma <- soma[meta$Treatment %in% c("Chow Vehicle", "DIO-NASH Vehicle"),]
# core_groups <- comb_outrm_meta$Treatment %in% c("Chow Vehicle", "DIO-NASH Vehicle")
# PC <- prcomp(comb_outrm2[core_groups,], center = TRUE, scale = TRUE)
# PC4 <- PC$rotation[,4]
# top_PC4_proteins <- names(PC4[order(abs(PC4), decreasing = TRUE)][1:100])
# 
# # Let's see if they're significantly regulated in one, the other, or in both datasets
# source("~/JDJG/scripts/signif_targets_function.R")
# JYNR210602_sig <- get_signif_targets(soma1 = outrm2_JYNR210602, meta1 = outrm_meta_JYNR210602,
#                                      group1 = 3, group2 = 1)
# JYNR210602_sig <- JYNR210602_sig[order(abs(JYNR210602_sig$log2fc), decreasing = TRUE),]
# JYNR210602_sig <- JYNR210602_sig[JYNR210602_sig$pval_adj < 0.05,]$seqid # 494 targets
# 
# 
# JYNR220601_sig <- get_signif_targets(soma1 = outrm2_JYNR220601, meta1 = outrm_meta_JYNR220601,
#                                      group1 = 2, group2 = 1)
# JYNR220601_sig <- JYNR220601_sig[order(abs(JYNR220601_sig$log2fc), decreasing = TRUE),]
# JYNR220601_sig <- JYNR220601_sig[JYNR220601_sig$pval_adj < 0.05,]$seqid # 3622 targets
# 
# ## SANITY CHECK
# # Top 100 PC4 proteins: Are they significant?
# sum(top_PC4_proteins[1:50] %in% JYNR210602_sig) # 9 (top50), 17 (top100)
# sum(top_PC4_proteins[1:50] %in% JYNR220601_sig) # 33 (top50), 67 (top100)
# 
# # Are they in the Top 100 proteins by FC?
# sum(top_PC4_proteins[1:50] %in% JYNR210602_sig[1:100]) # 0 (top50), 1 (top100)
# sum(top_PC4_proteins[1:50] %in% JYNR220601_sig[1:100]) # 0 (top50), 1 (top100)
# # No.
# 
# # Are they in the Top 200 proteins by FC?
# sum(top_PC4_proteins %in% JYNR210602_sig[1:200]) # 4
# sum(top_PC4_proteins %in% JYNR220601_sig[1:200]) # 6
# # Not really.
# 
# # Plot it (volcano style)
# myplots <- list(volcano_plot(soma1 = outrm2_JYNR210602, meta1 = outrm_meta_JYNR210602,
#                              custom_col_targets = top_PC4_proteins[1:50], custom_col_name = "PC4 Top50",
#                              title = "JYNR210602, case vs. control"),
#                 volcano_plot(soma1 = outrm2_JYNR220601, meta1 = outrm_meta_JYNR220601, group1 = 2,
#                              custom_col_targets = top_PC4_proteins[1:50], custom_col_name = "PC4 Top50",
#                              title = "JYNR220601, case vs. control"),
#                 volcano_plot(soma1 = outrm2_JYNR210602, meta1 = outrm_meta_JYNR210602,
#                              soma2 = outrm_JYNR220601, meta2 = outrm_meta_JYNR220601,
#                              group1 = 1, group2 = 1,
#                              custom_col_targets = top_PC4_proteins[1:50], custom_col_name = "PC4 Top50",
#                              title = "Both studies, Chow Vehicle"),
#                 volcano_plot(soma1 = outrm2_JYNR210602, meta1 = outrm_meta_JYNR210602,
#                              soma2 = outrm_JYNR220601, meta2 = outrm_meta_JYNR220601,
#                              group1 = 3, group2 = 2,
#                              custom_col_targets = top_PC4_proteins[1:50], custom_col_name = "PC4 Top50",
#                              title = "Both studies, DIO-NASH Vehicle"))
# legend <- get_legend(volcano_plot(legend.dir = "horizontal",
#                                   custom_col_targets = top_PC4_proteins[1:50], custom_col_name = "PC4 Top50"))
# g <- ggarrange(plotlist = myplots,
#                ncol = 2, nrow = 2,
#                common.legend = TRUE,
#                legend.grob = legend,
#                hjust = -2)
# title <- "Volcano plots, outliers removed, PC4 top 50 proteins marked"
# filename <- "figures/JYNR210602_JYNR220601_volc_PC4_top50"
# annotate_figure(g, top = text_grob(title))
# filename <- paste0(filename, ".png")
# ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 185, bg = "white")
# # This is just a sanity check. The proteins driving the inter-study differences would
# # be expected to be highly significant when contrasting between studies.
# # It does not seem to be the case (at least not as convincingly as I would have expected).
# # Edit: For top 50 it is actually alright. For intra-study analyses, they are rarely significant,
# # and for inter-study analyses they are often significant.
# 
# # Next step: Check if they're stable in calibrator samples.
# ## Create version with outliers removed prior to normalisation
# # Load normalisation script
# source("scripts/soma_norm_function.R")
# 
# # Load unnormalised adats
# adat1_unnorm <- read_adat("~/JDJG/data/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.adat")
# adat2_unnorm <- read_adat("~/JDJG/data/raw/Mouse_studies/JYNR220601_GUS2022_898/SS-2231349_v4.1_other.adat")
# adat1_norm1 <- mouse_JYNR210602 <- read_adat("~/JDJG/data/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.hybNorm.medNormInt.plateScale.adat")
# adat2_norm1 <- mouse_JYNR220601 <- read_adat("~/JDJG/data/raw/Mouse_studies/JYNR220601_GUS2022_898/SS-2231349_v4.1_other.hybNorm.medNormInt.plateScale.adat")
# 
# # Remove outliers
# out1 <- c(4018697918, 4018696639, 4018711275, 4018700427)
# keep1 <- !adat1_unnorm$SampleId %in% out1
# unnorm1 <- adat1_unnorm[keep1,]; anno1 <- adat1_norm1[keep1,]
# 
# out2 <- c(118606)
# keep2 <- !adat2_unnorm$SampleId %in% out2
# unnorm2 <- adat2_unnorm[keep2,]; anno2 <- adat2_norm1[keep2,]
# 
# # Normalisation (steps 1-4)
# temp_outrm2_JYNR210602 <- soma_norm(adat = unnorm1,
#                                     annotation_df = anno1,
#                                     medNormSMP = FALSE)
# 
# temp_outrm2_JYNR220601 <- soma_norm(adat = unnorm2,
#                                     annotation_df = anno2,
#                                     medNormSMP = FALSE)
# 
# # Select calibrator samples
# cal_JYNR210602 <- temp_outrm2_JYNR210602[temp_outrm2_JYNR210602$SampleType == "Calibrator",]
# cal_JYNR210602 <- cal_JYNR210602[, colnames(cal_JYNR210602) %in% colnames(outrm2_JYNR210602)]
# cal_JYNR210602 <- log2(cal_JYNR210602)
# 
# cal_JYNR220601 <- temp_outrm2_JYNR220601[temp_outrm2_JYNR220601$SampleType == "Calibrator",]
# cal_JYNR220601 <- cal_JYNR220601[, colnames(cal_JYNR220601) %in% colnames(outrm2_JYNR220601)]
# cal_JYNR220601 <- log2(cal_JYNR220601)
# 
# 
# # Unnormalised calibrator samples selection
# cal_unnorm_JYNR210602 <- adat1_unnorm[adat1_unnorm$SampleType == "Calibrator",]
# cal_unnorm_JYNR210602 <- cal_unnorm_JYNR210602[, colnames(cal_unnorm_JYNR210602) %in% colnames(outrm2_JYNR210602)]
# cal_unnorm_JYNR210602 <- log2(cal_unnorm_JYNR210602)
# 
# cal_unnorm_JYNR220601 <- adat2_unnorm[adat2_unnorm$SampleType == "Calibrator",]
# cal_unnorm_JYNR220601 <- cal_unnorm_JYNR220601[, colnames(cal_unnorm_JYNR220601) %in% colnames(outrm2_JYNR220601)]
# cal_unnorm_JYNR220601 <- log2(cal_unnorm_JYNR220601)
# 
# # Calibrator volcano plot
# myplots <- list(volcano_plot(soma1 = cal_unnorm_JYNR210602, soma2 = cal_unnorm_JYNR220601, meta_aligned = TRUE,
#                              title = "Unnormalised",
#                              custom_col_targets = top_PC4_proteins[1:50], custom_col_name = "PC4 Top50"),
#                 volcano_plot(soma1 = cal_JYNR210602, soma2 = cal_JYNR220601, meta_aligned = TRUE,
#                              title = "Normalised",
#                              custom_col_targets = top_PC4_proteins[1:50], custom_col_name = "PC4 Top50"))
# g <- ggarrange(plotlist = myplots,
#                ncol = 2, nrow = 1,
#                common.legend = TRUE,
#                legend.grob = legend,
#                hjust = -2)
# title <- "Calibrator samples, JYNR210602 vs. JYNR220601, PC4 top 50 proteins marked"
# filename <- "figures/JYNR210602_JYNR220601_volc_PC4_top50_cal"
# annotate_figure(g, top = text_grob(title))
# filename <- paste0(filename, ".png")
# ggsave(filename, width = 4000, height = 1200, units = "px", dpi = 185, bg = "white")
# 
# # Calibrator volcano plot (top 10)
# myplots <- list(volcano_plot(soma1 = cal_unnorm_JYNR210602, soma2 = cal_unnorm_JYNR220601, meta_aligned = TRUE,
#                              title = "Unnormalised",
#                              custom_col_targets = top_PC4_proteins[1:10], custom_col_name = "PC4 Top10"),
#                 volcano_plot(soma1 = cal_JYNR210602, soma2 = cal_JYNR220601, meta_aligned = TRUE,
#                              title = "Normalised",
#                              custom_col_targets = top_PC4_proteins[1:10], custom_col_name = "PC4 Top10"))
# g <- ggarrange(plotlist = myplots,
#                ncol = 2, nrow = 1,
#                common.legend = TRUE,
#                legend.grob = legend,
#                hjust = -2)
# title <- "Calibrator samples, JYNR210602 vs. JYNR220601, PC4 top 10 proteins marked"
# filename <- "figures/JYNR210602_JYNR220601_volc_PC4_top10_cal"
# annotate_figure(g, top = text_grob(title))
# filename <- paste0(filename, ".png")
# ggsave(filename, width = 4000, height = 1200, units = "px", dpi = 185, bg = "white")
# 
# 
# # Data standardisation exploration (volcano plots)
# # Define functions
# cal_norm <- function(adat = unnorm1, annotation_df = anno1,
#                      hybNorm = TRUE, medNormInt = TRUE,
#                      plateScale = TRUE, medNormSMP = FALSE){
#   cal_df <- soma_norm(adat = adat, annotation_df = annotation_df,
#                       hybNorm = hybNorm, medNormInt = medNormInt,
#                       plateScale = plateScale, medNormSMP = medNormSMP)
#   cal_df <- cal_df[cal_df$SampleType == "Calibrator",]
#   cal_df <- cal_df[, colnames(cal_df) %in% colnames(outrm2_JYNR210602)]
#   cal_df <- log2(cal_df)
#   return (cal_df)
# }
# cal_norm_both <- function(hybNorm = FALSE, medNormInt = FALSE,
#                           plateScale = FALSE, medNormSMP = FALSE){
#   cal_df1 <- cal_norm(hybNorm = hybNorm, medNormInt = medNormInt,
#                       plateScale = plateScale, medNormSMP = medNormSMP)
#   cal_df2 <- cal_norm(adat = unnorm2, annotation_df = anno2,
#                       hybNorm = hybNorm, medNormInt = medNormInt,
#                       plateScale = plateScale, medNormSMP = medNormSMP)
#   cal_dfs <- list(cal_df1, cal_df2)
#   return (cal_dfs)
# }
# 
# cal_hybNorm <- cal_norm_both(hybNorm = TRUE)
# cal_medInt <- cal_norm_both(hybNorm = TRUE, medNormInt = TRUE)
# cal_plateScale <- cal_norm_both(hybNorm = TRUE, medNormInt = TRUE,
#                                 plateScale = TRUE)
# # Du fik en fejl her. Gå op og kør normaliseringsfunktionerne separat og find ud af, hvad problemet er.
# # Hvis du ikke når det inden mødet, kan du nøjes med hurtigt at smide en medNormInt volcano ind.
# cal_medSMP <- cal_norm_both(hybNorm = TRUE, medNormInt = TRUE,
#                             plateScale = TRUE, medNormSMP = TRUE)
# 
# myplots <- list(volcano_plot(soma1 = cal_unnorm_JYNR210602, soma2 = cal_unnorm_JYNR220601, meta_aligned = TRUE,
#                              title = "Unnormalised",
#                              custom_col_targets = top_PC4_proteins[1:50], custom_col_name = "PC4 Top50"),
#                 volcano_plot(soma1 = cal_JYNR210602, soma2 = cal_JYNR220601, meta_aligned = TRUE,
#                              title = "Normalised",
#                              custom_col_targets = top_PC4_proteins[1:50], custom_col_name = "PC4 Top50"))
# g <- ggarrange(plotlist = myplots,
#                ncol = 2, nrow = 1,
#                common.legend = TRUE,
#                legend.grob = legend,
#                hjust = -2)
# title <- "Calibrator samples, JYNR210602 vs. JYNR220601, PC4 top 50 proteins marked"
# filename <- "figures/JYNR210602_JYNR220601_volc_PC4_top50_cal"
# annotate_figure(g, top = text_grob(title))
# filename <- paste0(filename, ".png")
# ggsave(filename, width = 4000, height = 1200, units = "px", dpi = 185, bg = "white")
# 
# cal_test_JYNR210602 <- cal_norm(adat = unnorm1, annotation_df = anno1,
#                                 hybNorm = TRUE, medNormInt = TRUE,
#                                 plateScale = FALSE, medNormSMP = FALSE)
# volcano_plot(soma1 = cal_dfs[[1]], soma2 = cal_dfs[[2]], meta_aligned = TRUE,
#              title = "Unnormalised",
#              custom_col_targets = top_PC4_proteins[1:10], custom_col_name = "PC4 Top10")
# 
# cal_test_JYNR220601 <- soma_norm(adat = unnorm2, annotation_df = anno2,
#                                  hybNorm = TRUE, medNormInt = TRUE,
#                                  plateScale = FALSE, medNormSMP = FALSE)
# 
# cal_test_JYNR210602 <- cal_test_JYNR210602[cal_test_JYNR210602$SampleType == "Calibrator",]
# cal_test_JYNR210602 <- cal_test_JYNR210602[, colnames(cal_test_JYNR210602) %in% colnames(outrm2_JYNR210602)]
# cal_test_JYNR210602 <- log2(cal_test_JYNR210602)
# 
# cal_test_JYNR220601 <- cal_test_JYNR220601[cal_test_JYNR220601$SampleType == "Calibrator",]
# cal_test_JYNR220601 <- cal_test_JYNR220601[, colnames(cal_test_JYNR220601) %in% colnames(outrm2_JYNR220601)]
# cal_test_JYNR220601 <- log2(cal_test_JYNR220601)
# 
# volcano_plot(soma1 = cal_test_JYNR210602, soma2 = cal_test_JYNR220601, meta_aligned = TRUE,
#              title = "Unnormalised",
#              custom_col_targets = top_PC4_proteins[1:10], custom_col_name = "PC4 Top10")
# 
# 
# # Calibrator density plot
# cals_unnorm <- merge_dfs(cal_unnorm_JYNR210602, cal_unnorm_JYNR220601)
# cals <- merge_dfs(cal_JYNR210602, cal_JYNR220601)
# cals$SampleId <- cals_unnorm$SampleId <- rownames(cals)
# cals$Dataset <- cals_unnorm$Dataset <- c(rep("JYNR210602",3), rep("JYNR220601", 3))
# cals_unnorm_long <- reshape2::melt(cals_unnorm, id.vars = c("SampleId", "Dataset"))
# cals_long <- reshape2::melt(cals, id.vars = c("SampleId", "Dataset"))
# cols <- rep("#619CFF", length(cals_long$Dataset))
# cols[cals_long$Dataset == "JYNR210602"] <- "#F8766D"
# png(filename = "~/JDJG/figures/adat_dens_cal.png", width = 2000, height = 400, units = "px", pointsize = 20)
# layout(mat = matrix(c(1,2), nrow = 1, ncol = 2, byrow=TRUE))
# density_compare_function(cals_unnorm_long$value, cals_unnorm_long$SampleId,
#                          lwd = 4, lty = 1, col = cols, alpha = 0.7, leg = unique(cals_unnorm_long$Dataset),
#                          xlab = "Log2-transformed aptamer measurements (RFU)")
# title(main = "Density distributions for calibrator samples (unnormalised)")
# density_compare_function(cals_long$value, cals_long$SampleId,
#                          lwd = 4, lty = 1, col = cols, alpha = 0.7, leg = unique(cals_long$Dataset),
#                          xlab = "Log2-transformed aptamer measurements (RFU)")
# title(main = "Density distributions for calibrator samples (normalised)")
# dev.off()
# 
# 
# # Calibrator density plot (PC4 top50 proteins)
# cals_top100 <- merge_dfs(cal_JYNR210602, cal_JYNR220601)
# cals_top100 <- cals_top100[,colnames(cals_top100) %in% top_PC4_proteins[1:100]]
# cals_top100$SampleId <- cals_unnorm$SampleId <- rownames(cals)
# cals_top100$Dataset <- cals_unnorm$Dataset <- c(rep("JYNR210602",3), rep("JYNR220601", 3))
# cals_top100_long <- reshape2::melt(cals_top100, id.vars = c("SampleId", "Dataset"))
# cols <- rep("#619CFF", length(cals_long$Dataset))
# cols[cals_top100_long$Dataset == "JYNR210602"] <- "#F8766D"
# cals_top50_long <- cals_top100_long[cals_top100_long$variable %in% top_PC4_proteins[1:50],]
# cals_top10_long <- cals_top100_long[cals_top100_long$variable %in% top_PC4_proteins[1:10],]
# png(filename = "~/JDJG/figures/adat_dens_cal_PC4top.png", width = 2000, height = 400, units = "px", pointsize = 20)
# layout(mat = matrix(c(1,2,3), nrow = 1, ncol = 3, byrow=TRUE))
# density_compare_function(cals_top100_long$value, cals_top100_long$SampleId,
#                          lwd = 4, lty = 1, col = cols, alpha = 0.7, leg = unique(cals_top100_long$Dataset),
#                          xlab = "Log2-transformed aptamer measurements (RFU)")
# title(main = "Calibrator samples, PC4 top 100 proteins")
# density_compare_function(cals_top50_long$value, cals_top50_long$SampleId,
#                          lwd = 4, lty = 1, col = cols, alpha = 0.7, leg = unique(cals_top50_long$Dataset),
#                          xlab = "Log2-transformed aptamer measurements (RFU)")
# title(main = "Calibrator samples, PC4 top 50 proteins")
# density_compare_function(cals_top10_long$value, cals_top10_long$SampleId,
#                          lwd = 4, lty = 1, col = cols, alpha = 0.7, leg = unique(cals_top10_long$Dataset),
#                          xlab = "Log2-transformed aptamer measurements (RFU)")
# title(main = "Calibrator samples, PC4 top 10 proteins")
# dev.off()
# 
# # Sanity check (randomly sample 10, 50 and 100 proteins to see what the distributions might look like)
# set.seed(421)
# rand_10 <- sample(1:dim(cal_JYNR210602)[2], size = 10, replace = FALSE)
# rand_50 <- sample(1:dim(cal_JYNR210602)[2], size = 50, replace = FALSE)
# rand_100 <- sample(1:dim(cal_JYNR210602)[2], size = 100, replace = FALSE)
# cals_rand10_long <- cals_long[cals_long$variable %in% colnames(cal_JYNR210602)[rand_10],]
# cals_rand50_long <- cals_long[cals_long$variable %in% colnames(cal_JYNR210602)[rand_50],]
# cals_rand100_long <- cals_long[cals_long$variable %in% colnames(cal_JYNR210602)[rand_100],]
# png(filename = "~/JDJG/figures/adat_dens_cal_rand.png", width = 2000, height = 400, units = "px", pointsize = 20)
# layout(mat = matrix(c(1,2,3), nrow = 1, ncol = 3, byrow=TRUE))
# density_compare_function(cals_rand100_long$value, cals_rand100_long$SampleId,
#                          lwd = 4, lty = 1, col = cols, alpha = 0.7, leg = unique(cals_rand100_long$Dataset),
#                          xlab = "Log2-transformed aptamer measurements (RFU)")
# title(main = "Calibrator samples, 100 random proteins")
# density_compare_function(cals_rand50_long$value, cals_rand50_long$SampleId,
#                          lwd = 4, lty = 1, col = cols, alpha = 0.7, leg = unique(cals_rand50_long$Dataset),
#                          xlab = "Log2-transformed aptamer measurements (RFU)")
# title(main = "Calibrator samples, 50 random proteins")
# density_compare_function(cals_rand10_long$value, cals_rand10_long$SampleId,
#                          lwd = 4, lty = 1, col = cols, alpha = 0.7, leg = unique(cals_rand10_long$Dataset),
#                          xlab = "Log2-transformed aptamer measurements (RFU)")
# title(main = "Calibrator samples, 10 random proteins")
# dev.off()
# cals_long
# 
# # from wrangle script:
# cal_JYNR210602 <- log2(cal_JYNR210602); cal_JYNR220601 <- log2(cal_JYNR220601)
# 
# 
# 
# 
# 
# 
# 
# ################################################################################
# ## Batch effect plot
# # Regular plot (baseline, limma, 1-5, cross-norm)
# png("figures/batch_effect_log.png", width = 2000, height = 1200,
#     pointsize = 28)
# col <- c(rep("#F8766D",dim(log_JYNR210602)[1]),rep("#00BFC4",dim(log_JYNR220601)[1]))
# batch <- comb_meta$dataset
# m <- matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
# layout(mat = m, heights = c(1,1,0.104)); par(mar = c(2,2,2,1))
# # m <- matrix(c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE)
# # layout(mat = m, heights = c(1,0.08)); par(mar = c(2,2,2,1))
# 
# 
# # Baseline
# boxplot(t(comb_log),
#         col=col,
#         main="Baseline")
# boxplot(t(comb_limma),
#         col=col,
#         main="Batch effects corrected (limma)")
# boxplot(t(comb_log_1_5),
#         col=col,
#         main="Median signal normalisation (SomaLogic)")
# boxplot(t(comb_crossnorm),
#         col=col,
#         main="Cross-normalised (SomaLogic)")
# 
# par(mar = c(0,0,0.5,0))
# plot(1, type="n", axes=FALSE, xlab="", ylab="")
# legend(x = "top",
#        legend = c("JYNR210602", "JYNR220601"),
#        fill = c("#F8766D", "#00BFC4"),
#        horiz=TRUE)
# dev.off()
# 
# png("figures/batch_effect_log.png", width = 2000, height = 1200,
#     pointsize = 28)
# col <- c(rep("#F8766D",dim(log_JYNR210602)[1]),rep("#00BFC4",dim(log_JYNR220601)[1]))
# batch <- comb_meta$dataset
# m <- matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
# layout(mat = m, heights = c(1,1,0.104)); par(mar = c(2,2,2,1))
# # m <- matrix(c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE)
# # layout(mat = m, heights = c(1,0.08)); par(mar = c(2,2,2,1))
# 
# 
# # Z-scores
# png("figures/batch_effect_z.png", width = 2000, height = 1200,
#     pointsize = 28)
# col <- c(rep("#F8766D",dim(log_JYNR210602)[1]),rep("#00BFC4",dim(log_JYNR220601)[1]))
# batch <- comb_meta$dataset
# m <- matrix(c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE)
# layout(mat = m, heights = c(1,0.08)); par(mar = c(2,2,2,1))
# boxplot(t(comb_z),
#         col=col,
#         main="Z-scores")
# boxplot(t(comb_z_log),
#         col=col,
#         main="Z-scores (from log2 data)")
# 
# par(mar = c(0,0,0.5,0))
# plot(1, type="n", axes=FALSE, xlab="", ylab="")
# legend(x = "top",
#        legend = c("JYNR210602", "JYNR220601"),
#        fill = c("#F8766D", "#00BFC4"),
#        horiz=TRUE)
# dev.off()
