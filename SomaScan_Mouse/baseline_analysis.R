## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel)
setwd("~/JDJG/clean/")

## Load data and perform initial wrangling (this should be further automated for future use)
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")
source("scripts/SomaScan_Mouse/CV_plot_function.R")
source("scripts/SomaScan_Mouse/overlap_function.R")


# ################################################################################
# ### Analysis
# ## Top 100 genes
top100_overlap <- function(soma1, soma2, meta1, meta2){
  ## Dataset1
  # Compute significance and fold change
  sig_tar_1 <- get_signif_targets(soma1 = soma1, meta1 = meta1, group1 = 3, group2 = 1)
  # Order by log2(FC) (absolute)
  sig_tar_1 <- sig_tar_1[order(abs(sig_tar_1$log2fc), decreasing = TRUE),]
  # Extract seqids for only significant targets
  sig_tar_1 <- sig_tar_1[sig_tar_1$pval_adj < 0.05,]$seqid
  
  ## Dataset2
  # Compute significance and fold change
  sig_tar_2 <- get_signif_targets(soma1 = soma2, meta1 = meta2, group1 = 2, group2 = 1)
  # Order by log2(FC) (absolute)
  sig_tar_2 <- sig_tar_2[order(abs(sig_tar_2$log2fc), decreasing = TRUE),]
  # Extract seqids for only significant targets
  sig_tar_2 <- sig_tar_2[sig_tar_2$pval_adj < 0.05,]$seqid
  
  # Get top100 overlap
  overlap <- sum(sig_tar_1[1:100] %in% sig_tar_2[1:100])
  
  return(overlap)
}



# JYNR210602_sig <- get_signif_targets(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
#                                      group1 = 3, group2 = 1)
# JYNR210602_sig$log2fc <- abs(JYNR210602_sig$log2fc)
# JYNR210602_sig <- JYNR210602_sig[order(abs(JYNR210602_sig$log2fc), decreasing = TRUE),]
# top_JYNR210602 <- JYNR210602_sig[JYNR210602_sig$pval_adj < 0.05,]$seqid
# top_JYNR210602 <- JYNR210602_sig$seqid
# 
# JYNR220601_sig <- get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
#                                      group1 = 2, group2 = 1)
# JYNR220601_sig$log2fc <- abs(JYNR220601_sig$log2fc)
# JYNR220601_sig <- JYNR220601_sig[order(JYNR220601_sig$log2fc, decreasing = TRUE),]
# top_JYNR220601 <- JYNR220601_sig[JYNR220601_sig$pval_adj < 0.05,]$seqid
# top_JYNR220601 <- JYNR220601_sig$seqid
# 
# # Top 10 in Top 10
# sum(top_JYNR210602[1:10] %in% top_JYNR220601[1:10]) # 55 (per seqid), w/o sig. threshold: 59 (per gene), 58 (per seqid)
# 
# # Top 100 in Top 100
# sum(top_JYNR210602[1:100] %in% top_JYNR220601[1:100]) # 55 (per seqid), w/o sig. threshold: 59 (per gene), 58 (per seqid)
# 
# # Top 100 in Top 200
# sum(top_JYNR210602[1:100] %in% top_JYNR220601[1:200]) # 84
# sum(top_JYNR220601[1:100] %in% top_JYNR210602[1:200]) # 59
# 
# # Top 100 in Top 500
# sum(top_JYNR210602[1:100] %in% top_JYNR220601[1:500]) # 99
# sum(top_JYNR220601[1:100] %in% top_JYNR210602[1:500]) # 59
# 
# # Top 100 in Top 1000
# sum(top_JYNR210602[1:100] %in% top_JYNR220601[1:1000]) # 99 (per seqid); w/o sig. 86 (per gene), 85 (per seqid)
# sum(top_JYNR220601[1:100] %in% top_JYNR210602[1:1000]) # 59 (clearly, 41/100 are not significant in JYNR210602)
# 
# 
# 
# ## Top 100 genes (with outliers removed)
# JYNR210602_sig <- get_signif_targets(soma1 = outrm2_JYNR210602, meta1 = outrm_meta_JYNR210602,
#                                      group1 = 3, group2 = 1)
# JYNR210602_sig$log2fc <- abs(JYNR210602_sig$log2fc)
# JYNR210602_sig <- JYNR210602_sig[order(abs(JYNR210602_sig$log2fc), decreasing = TRUE),]
# top_JYNR210602 <- JYNR210602_sig[JYNR210602_sig$pval_adj < 0.05,]$seqid
# top_JYNR210602 <- JYNR210602_sig$seqid
# 
# JYNR220601_sig <- get_signif_targets(soma1 = outrm2_JYNR220601, meta1 = outrm_meta_JYNR220601,
#                                      group1 = 2, group2 = 1)
# JYNR220601_sig$log2fc <- abs(JYNR220601_sig$log2fc)
# JYNR220601_sig <- JYNR220601_sig[order(JYNR220601_sig$log2fc, decreasing = TRUE),]
# top_JYNR220601 <- JYNR220601_sig[JYNR220601_sig$pval_adj < 0.05,]$seqid
# top_JYNR220601 <- JYNR220601_sig$seqid
# 
# # Top 10 in Top 10
# sum(top_JYNR210602[1:10] %in% top_JYNR220601[1:10]) # 9
# 
# # Top 100 in Top 100
# sum(top_JYNR210602[1:100] %in% top_JYNR220601[1:100]) # 75
# 
# # Top 100 in Top 200
# sum(top_JYNR210602[1:100] %in% top_JYNR220601[1:200]) # 93
# sum(top_JYNR220601[1:100] %in% top_JYNR210602[1:200]) # 84
# 
# # Top 100 in Top 500
# sum(top_JYNR210602[1:100] %in% top_JYNR220601[1:500]) # 99
# sum(top_JYNR220601[1:100] %in% top_JYNR210602[1:500]) # 84
# 
# # Top 100 in Top 1000
# sum(top_JYNR210602[1:100] %in% top_JYNR220601[1:1000]) # 99
# sum(top_JYNR220601[1:100] %in% top_JYNR210602[1:1000]) # 84 (clearly, 16/100 are not significant in JYNR210602)



################################################################################
#### CV plots!


## CV% for different dfs
# # Baseline and outlier removal
# cv_plot(init_JYNR210602, init_JYNR220601, init_meta_JYNR210602, init_meta_JYNR220601,
#         text = "no outliers removed")
# ggsave("figures/cv_init.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# cv_plot(old_log_JYNR210602, old_log_JYNR220601, old_meta_JYNR210602, old_meta_JYNR220601,
#         text = "two outliers removed")
# ggsave("figures/cv_old.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# cv_plot(log_JYNR210602, log_JYNR220601, meta_mouse_JYNR210602, meta_mouse_JYNR220601,
#         text = "outliers removed")
# ggsave("figures/cv_base.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# 
# # Baseline and outlier removal (steps 1-5)
# cv_plot(init_1_5_JYNR210602, init_1_5_JYNR220601, init_meta_JYNR210602, init_meta_JYNR220601,
#         text = "no outliers removed, interplate median signal normalisation")
# ggsave("figures/cv_1_5_init.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# cv_plot(old_1_5_JYNR210602, old_1_5_JYNR220601, old_meta_JYNR210602, old_meta_JYNR220601,
#         text = "two outliers removed, interplate median signal normalisation")
# ggsave("figures/cv_1_5_old.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
# cv_plot(log_1_5_JYNR210602, log_1_5_JYNR220601, meta_mouse_JYNR210602, meta_mouse_JYNR220601,
#         text = "outliers removed, interplate median signal normalisation")
# ggsave("figures/cv_1_5_base.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")

# Thesis version
myplots <- list(cv_plot(init_JYNR210602, init_JYNR220601, init_meta_JYNR210602, init_meta_JYNR220601,
                        text = "No outliers removed", only_cont = TRUE, x_max = 15,
                        custom_title = TRUE, fig1 = TRUE),
                cv_plot(init_1_5_JYNR210602, init_1_5_JYNR220601, init_meta_JYNR210602, init_meta_JYNR220601,
                        text = "No outliers removed, interplate median signal normalisation",
                        only_cont = TRUE, x_max = 15, custom_title = TRUE, fig1 = TRUE),
                cv_plot(log_JYNR210602, log_JYNR220601, meta_mouse_JYNR210602, meta_mouse_JYNR220601,
                        text = "Outliers removed", only_cont = TRUE, x_max = 15,
                        custom_title = TRUE, fig1 = TRUE)
)
legend <- cv_plot(init_JYNR210602, init_JYNR220601, init_meta_JYNR210602, init_meta_JYNR220601,
                  get_legend = TRUE, only_cont = TRUE)
cv_final <- ggarrange(plotlist = myplots,
                      nrow = length(myplots),
                      common.legend = TRUE,
                      legend.grob = legend,
                      legend = "right",
                      hjust = -2)
# cv_final
# ggsave("figures/cv_final.png", width = 1000, height = 2000, units = "px", dpi = 150, bg = "white")
# 4000 x 2000
# halvdelen af det. Halvdelen af bredden.
# 2000 x 2000.
# så kvadratisk.


# Skip the rest of these CV plots




# ## Coefficient of variation (CV%)
# cvs <- matrix(nrow = dim(log_JYNR210602)[2],
#               ncol = 4)
# cvs[,1] <- sapply(1:ncol(log_JYNR210602), function(x) sd(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,x]) / mean(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,x]) * 100)
# cvs[,2] <- sapply(1:ncol(log_JYNR220601), function(x) sd(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,x]) / mean(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,x]) * 100)
# cvs[,3] <- sapply(1:ncol(log_JYNR210602), function(x) sd(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,x]) / mean(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,x]) * 100)
# cvs[,4] <- sapply(1:ncol(log_JYNR220601), function(x) sd(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,x]) / mean(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,x]) * 100)
# cvs <- as.data.frame(cvs); names(cvs) <- c("JYNR210602_placebo", "JYNR220601_placebo", "JYNR210602_case", "JYNR220601_case")
# 
# xmin <- min(cvs); xmax <- max(cvs)
# 
# p1 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_placebo, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_placebo, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Coefficient of variation (%) per aptamer") +
#   ggtitle("Placebo vs. placebo") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#                             legend.direction = "horizontal")
# 
# p2 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_case, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_case, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Coefficient of variation (%) per aptamer") +
#   ggtitle("Case vs. case") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# legend <- get_legend(p1)
# ggarrange(p1,p2,
#           common.legend = TRUE,
#           legend.grob = legend,
#           hjust = -2)
# ggsave("figures/mouse_cv.png", width = 2000, height = 1000, units = "px", dpi = 200, bg = "white")
# 
# 
# 
# ## Coefficient of variation (CV%) for the step 1-5 dataset
# cvs <- matrix(nrow = dim(log_1_5_JYNR210602)[2],
#               ncol = 4)
# cvs[,1] <- sapply(1:ncol(log_1_5_JYNR210602), function(x) sd(log_1_5_JYNR210602[meta_mouse_JYNR210602$Group == 1,x]) / mean(log_1_5_JYNR210602[meta_mouse_JYNR210602$Group == 1,x]) * 100)
# cvs[,2] <- sapply(1:ncol(log_1_5_JYNR220601), function(x) sd(log_1_5_JYNR220601[meta_mouse_JYNR220601$Group == 1,x]) / mean(log_1_5_JYNR220601[meta_mouse_JYNR220601$Group == 1,x]) * 100)
# cvs[,3] <- sapply(1:ncol(log_1_5_JYNR210602), function(x) sd(log_1_5_JYNR210602[meta_mouse_JYNR210602$Group == 3,x]) / mean(log_1_5_JYNR210602[meta_mouse_JYNR210602$Group == 3,x]) * 100)
# cvs[,4] <- sapply(1:ncol(log_1_5_JYNR220601), function(x) sd(log_1_5_JYNR220601[meta_mouse_JYNR220601$Group == 2,x]) / mean(log_1_5_JYNR220601[meta_mouse_JYNR220601$Group == 2,x]) * 100)
# cvs <- as.data.frame(cvs); names(cvs) <- c("JYNR210602_placebo", "JYNR220601_placebo", "JYNR210602_case", "JYNR220601_case")
# 
# xmin <- min(cvs); xmax <- max(cvs)
# 
# p1 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_placebo, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_placebo, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Coefficient of variation (%) per aptamer") +
#   ggtitle("Placebo vs. placebo, steps 1-5") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# p2 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_case, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_case, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Coefficient of variation (%) per aptamer") +
#   ggtitle("Case vs. case, steps 1-5") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# legend <- get_legend(p1)
# ggarrange(p1,p2,
#           common.legend = TRUE,
#           legend.grob = legend,
#           hjust = -2)
# ggsave("figures/mouse_1_5_cv.png", width = 2000, height = 1000, units = "px", dpi = 200, bg = "white")
# 
# 
# 
# 
# ## Coefficient of variation (CV%) for the outlier removal data
# cvs <- matrix(nrow = dim(outrm_JYNR210602)[2],
#               ncol = 4)
# cvs[,1] <- sapply(1:ncol(log_1_5_JYNR210602), function(x) sd(outrm_JYNR210602[outrm_meta_JYNR210602$Group == 1,x]) / mean(outrm_JYNR210602[outrm_meta_JYNR210602$Group == 1,x]) * 100)
# cvs[,2] <- sapply(1:ncol(log_1_5_JYNR220601), function(x) sd(outrm_JYNR220601[outrm_meta_JYNR220601$Group == 1,x]) / mean(outrm_JYNR220601[outrm_meta_JYNR220601$Group == 1,x]) * 100)
# cvs[,3] <- sapply(1:ncol(log_1_5_JYNR210602), function(x) sd(outrm_JYNR210602[outrm_meta_JYNR210602$Group == 3,x]) / mean(outrm_JYNR210602[outrm_meta_JYNR210602$Group == 3,x]) * 100)
# cvs[,4] <- sapply(1:ncol(log_1_5_JYNR220601), function(x) sd(outrm_JYNR220601[outrm_meta_JYNR220601$Group == 2,x]) / mean(outrm_JYNR220601[outrm_meta_JYNR220601$Group == 2,x]) * 100)
# cvs <- as.data.frame(cvs); names(cvs) <- c("JYNR210602_placebo", "JYNR220601_placebo", "JYNR210602_case", "JYNR220601_case")
# 
# xmin <- min(cvs); xmax <- max(cvs)
# 
# p1 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_placebo, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_placebo, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Coefficient of variation (%) per aptamer") +
#   ggtitle("Placebo vs. placebo, outliers removed") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# p2 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_case, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_case, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Coefficient of variation (%) per aptamer") +
#   ggtitle("Case vs. case, outliers removed (post-norm)") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# legend <- get_legend(p1)
# ggarrange(p1,p2,
#           common.legend = TRUE,
#           legend.grob = legend,
#           hjust = -2)
# ggsave("figures/mouse_outrm_cv.png", width = 2000, height = 1000, units = "px", dpi = 200, bg = "white")
# 
# 
# 
# ## Coefficient of variation (CV%) for the outlier removal data, removed pre-normalisation
# cvs <- matrix(nrow = dim(outrm2_JYNR210602)[2],
#               ncol = 4)
# cvs[,1] <- sapply(1:ncol(log_1_5_JYNR210602), function(x) sd(outrm2_JYNR210602[outrm_meta_JYNR210602$Group == 1,x]) / mean(outrm2_JYNR210602[outrm_meta_JYNR210602$Group == 1,x]) * 100)
# cvs[,2] <- sapply(1:ncol(log_1_5_JYNR220601), function(x) sd(outrm2_JYNR220601[outrm_meta_JYNR220601$Group == 1,x]) / mean(outrm2_JYNR220601[outrm_meta_JYNR220601$Group == 1,x]) * 100)
# cvs[,3] <- sapply(1:ncol(log_1_5_JYNR210602), function(x) sd(outrm2_JYNR210602[outrm_meta_JYNR210602$Group == 3,x]) / mean(outrm2_JYNR210602[outrm_meta_JYNR210602$Group == 3,x]) * 100)
# cvs[,4] <- sapply(1:ncol(log_1_5_JYNR220601), function(x) sd(outrm2_JYNR220601[outrm_meta_JYNR220601$Group == 2,x]) / mean(outrm2_JYNR220601[outrm_meta_JYNR220601$Group == 2,x]) * 100)
# cvs <- as.data.frame(cvs); names(cvs) <- c("JYNR210602_placebo", "JYNR220601_placebo", "JYNR210602_case", "JYNR220601_case")
# 
# xmin <- min(cvs); xmax <- max(cvs)
# 
# p1 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_placebo, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_placebo, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Coefficient of variation (%) per aptamer") +
#   ggtitle("Placebo vs. placebo, outliers removed") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# p2 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_case, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_case, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Coefficient of variation (%) per aptamer") +
#   ggtitle("Case vs. case, outliers removed (pre-norm)") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# legend <- get_legend(p1)
# ggarrange(p1,p2,
#           common.legend = TRUE,
#           legend.grob = legend,
#           hjust = -2)
# ggsave("figures/mouse_outrm2_cv.png", width = 2000, height = 1000, units = "px", dpi = 200, bg = "white")
# 
# 
# 
# ## Coefficient of variation (CV%) for stepwise outlier removal
# cv_data <- data.frame(control = rep(NA,10),
#                       case = rep(NA, 10),
#                       dataset = c(rep("JYNR210602",9), "JYNR220601"),
#                       outlier_rm = c(0:8,0),
#                       case_grp = c(rep(3,9),2))
# 
# for (i in 1:dim(cv_data)[1]){
#   df <- get(paste0("log_", cv_data$dataset[i])); meta <- get(paste0("meta_mouse_", cv_data$dataset[i]))
#   case <- df[meta$Group == cv_data$case_grp[i],]; control <- df[meta$Group == 1,]
#   out <- cv_data$outlier_rm[i]
#   if (!out == 0){
#     case <- case[-out,]
#     if (out == 8){
#       control <- control[-c(2,6),]
#     }
#     else {
#       control <- control[-out,]
#     }
#   }
#   cv_data$case[i] <- median(sapply(1:ncol(df), function(x) sd(case[,x]) / mean(case[,x]) * 100))
#   cv_data$control[i] <- median(sapply(1:ncol(df), function(x) sd(control[,x]) / mean(control[,x]) * 100))
# }
# 
# 
# 
# ## MAD
# cvs <- matrix(nrow = dim(log_JYNR210602)[2],
#               ncol = 4)
# cvs[,1] <- sapply(1:ncol(log_JYNR210602), function(x) mad(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,x]))
# cvs[,2] <- sapply(1:ncol(log_JYNR220601), function(x) mad(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,x]))
# cvs[,3] <- sapply(1:ncol(log_JYNR210602), function(x) mad(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,x]))
# cvs[,4] <- sapply(1:ncol(log_JYNR220601), function(x) mad(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,x]))
# cvs <- as.data.frame(cvs); names(cvs) <- c("JYNR210602_placebo", "JYNR220601_placebo", "JYNR210602_case", "JYNR220601_case")
# 
# xmin <- min(cvs); xmax <- max(cvs)
# 
# p1 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_placebo, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_placebo, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Median absolute deviation per aptamer") +
#   ggtitle("Placebo vs. placebo") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# p2 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_case, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_case, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Median absolute deviation per aptamer") +
#   ggtitle("Case vs. case") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# legend <- get_legend(p1)
# ggarrange(p1,p2,
#           common.legend = TRUE,
#           legend.grob = legend,
#           hjust = -2)
# ggsave("figures/mouse_mad.png", width = 2000, height = 1000, units = "px", dpi = 200, bg = "white")
# 
# 
# 
# ## Relative coefficient of variation (%) per aptamer
# # (MAD * 1.4826 / median)
# # Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9196089/
# cvs <- matrix(nrow = dim(log_JYNR210602)[2],
#               ncol = 4)
# cvs[,1] <- sapply(1:ncol(log_JYNR210602), function(x) mad(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,x]) / median(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,x]) * 100)
# cvs[,2] <- sapply(1:ncol(log_JYNR220601), function(x) mad(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,x]) / median(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,x]) * 100)
# cvs[,3] <- sapply(1:ncol(log_JYNR210602), function(x) mad(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,x]) / median(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,x]) * 100)
# cvs[,4] <- sapply(1:ncol(log_JYNR220601), function(x) mad(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,x]) / median(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,x]) * 100)
# cvs <- as.data.frame(cvs); names(cvs) <- c("JYNR210602_placebo", "JYNR220601_placebo", "JYNR210602_case", "JYNR220601_case")
# 
# xmin <- min(cvs); xmax <- max(cvs)
# 
# p1 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_placebo, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_placebo, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Relative coefficient of variation (%) per aptamer") +
#   ggtitle("Placebo vs. placebo") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# p2 <- cvs %>% ggplot() +
#   geom_histogram(aes(x = JYNR210602_case, fill = "JYNR210602"), alpha = 0.8) +
#   geom_histogram(aes(x = JYNR220601_case, fill = "JYNR220601"), alpha = 0.8) +
#   xlab("Relative coefficient of variation (%) per aptamer") +
#   ggtitle("Case vs. case") +
#   xlim(c(xmin,xmax)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.direction = "horizontal")
# 
# legend <- get_legend(p1)
# ggarrange(p1,p2,
#           common.legend = TRUE,
#           legend.grob = legend,
#           hjust = -2)
# ggsave("figures/mouse_rcv_m.png", width = 2000, height = 1000, units = "px", dpi = 200, bg = "white")
