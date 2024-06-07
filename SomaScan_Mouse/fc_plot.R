## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel)
library(WGCNA); library(scales)
setwd("~/JDJG/clean")
source("scripts/SomaScan_Mouse/overlap_function.R")

## Load data and perform initial wrangling (this should be further automated for future use)
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")


################################################################################
### Analysis
## Define plotting function
fc_plot <- function(handle = "log_", text = "", width = 1400, height = 700,
                    top100 = TRUE, signif = TRUE, save = TRUE){
  
  # Set plot margins (reset because it is modified by other code snippets)
  par(mar = c(5, 4, 4, 2) + 0.1 - c(1,0,2,0))
  # low, left, top, right
  
  # Get data
  soma1 <- get(paste0(handle, "JYNR210602"))
  soma2 <- get(paste0(handle, "JYNR220601"))
  if ("log" %in% str_split(handle, "_")[[1]][1]){
    meta1 <- meta_mouse_JYNR210602
    meta2 <- meta_mouse_JYNR220601
  } else if ("init" %in% str_split(handle, "_")[[1]][1]){
    meta1 <- init_meta_JYNR210602
    meta2 <- init_meta_JYNR220601
  } else if ("old" %in% str_split(handle, "_")[[1]][1]){
    meta1 <- old_meta_JYNR210602
    meta2 <- old_meta_JYNR220601
  }
  
  ## Define plotting variables
  JYNR210602_sig <- get_signif_targets(soma1 = soma1, meta1 = meta1,
                                       group1 = 3, group2 = 1)
  
  JYNR220601_sig <- get_signif_targets(soma1 = soma2, meta1 = meta2,
                                       group1 = 2, group2 = 1)
  
  X <- JYNR210602_sig$log2fc
  Y <- JYNR220601_sig$log2fc
  
  cols <- c("#FFA355", "#55B1FF", "#78ED56")
  
  # Compute correlations and axis limits
  fit <- lm(Y ~ X)
  corr <- corAndPvalue(X,Y, method = "spearman")
  corr2 <- corAndPvalue(X,Y)
  xlim <- c(min(X)*1.1,max(X)*1.1); ylim <- c(min(Y)*1.1,max(Y)*1.1)
  
  # title <- paste0("FC vs. FC, ", text)
  title <- text
  
  ### Point plot of fold change (Soma) vs. fold change (MS)
  if (top100 == TRUE){
    # Colour by top 100 genes
    Top100_x <- rownames(JYNR210602_sig)[order(abs(X), decreasing = TRUE)][1:100]
    Top100_y <- rownames(JYNR220601_sig)[order(abs(Y), decreasing = TRUE)][1:100]
    Top100_col <- rep("lightgrey", length(X))
    Top100_col[rownames(JYNR210602_sig) %in% Top100_x] <- cols[1]
    Top100_col[rownames(JYNR210602_sig) %in% Top100_y] <- cols[2]
    Top100_col[rownames(JYNR210602_sig) %in% Top100_x & rownames(JYNR210602_sig) %in% Top100_y] <- cols[3]
    
    # Plot it
    if (save == TRUE){
      png(file=paste0("~/JDJG/figures/mouse_FC_", handle, "top100.png"),
          width = width, height = height, units = "px", pointsize = 28)
    }
    col <- Top100_col
    plot(X,
         Y,
         pch = 19, col = alpha(col, 0.8),
         xlim = xlim,
         ylim = ylim,
         main = title,
         xlab = "Fold change (SomaMouse1)", ylab = "Fold change (SomaMouse2)")
    lines(X,
          # X * coef(fit)[2] + coef(fit)[1],
          X,
          col = 1)
    rr <- vector("expression", 3)
    rr[1] <- substitute(expression(italic(R)^2 == val),
                        list(val = format(unlist(summary(fit)[9]), dig=2)))[2]
    rr[2] <- paste0("Pearson = ", format(corr2$cor, dig=2))
    rr[3] <- paste0("Spearman = ", format(corr$cor, dig=2))
    # rr <- rr[-1]
    legend("topleft", legend = rr, bty = "n")
    legend("bottomright", legend = c("Top100 in SomaMouse1",
                                     "Top100 in SomaMouse2",
                                     "Top100 in both"),
           col = cols, pch = 19)
    if (save == TRUE){
      dev.off()
    }
  }
  
  
  if (signif == TRUE){
    # Colour by significance (instead of top 100)
    # At the moment, the "sign" variable does not mean significance, but rather
    # significance + abs(log2fc) > 0.3. This is impractical, but I have too
    # much code depending on the get_signif_targets to want to change the way
    # it works. SO to get significance, don't select JYNR210602_sig$sign == TRUE,
    # but instead choose JYNR210602_sig$pval_adj < 0.05.
    sig_x <- rownames(JYNR210602_sig[JYNR210602_sig$pval_adj < 0.05,])
    sig_y <- rownames(JYNR220601_sig[JYNR220601_sig$pval_adj < 0.05,])
    sig_col <- rep("lightgrey", length(X))
    sig_col[rownames(JYNR210602_sig) %in% sig_x] <- cols[1]
    sig_col[rownames(JYNR210602_sig) %in% sig_y] <- cols[2]
    sig_col[rownames(JYNR210602_sig) %in% sig_x & rownames(JYNR210602_sig) %in% sig_y] <- cols[3]
    
    # Plot it
    if (save == TRUE){
      png(file=paste0("~/JDJG/figures/mouse_FC_", handle, "signif.png"),
          width = width, height = height, units = "px", pointsize = 28)
    }
    col <- sig_col
    plot(X,
         Y,
         pch = 19, col = alpha(col, 0.8),
         xlim = xlim,
         ylim = ylim,
         main = title,
         xlab = "Fold change (SomaMouse1)", ylab = "Fold change (SomaMouse2)")
    lines(X,
          # X * coef(fit)[2] + coef(fit)[1],
          X,
          col = 1)
    rr <- vector("expression", 3)
    rr[1] <- substitute(expression(italic(R)^2 == val),
                        list(val = format(unlist(summary(fit)[9]), dig=2)))[2]
    rr[2] <- paste0("Pearson = ", format(corr2$cor, dig=2))
    rr[3] <- paste0("Spearman = ", format(corr$cor, dig=2))
    # rr <- rr[-1]
    legend("topleft", legend = rr, bty = "n")
    legend("bottomright", legend = c("Significant in SomaMouse1",
                                     "Significant in SomaMouse2",
                                     "Significant in both"),
           col = cols, pch = 19)
    if (save == TRUE){
      dev.off()
    } 
  }
}

# ## Plotting
# # Baseline and outlier removal
# fc_plot(handle = "init_",
#         text = "no outliers removed")
# # fc_plot(handle = "old_log_",
# #         text = "two outliers removed")
# fc_plot(handle = "log_",
#         text = "outliers removed")
# 
# # Baseline and outlier removal (steps 1-5)
# fc_plot(handle = "init_1_5_",
#         text = "no outliers removed, medNorm")
# # fc_plot(handle = "old_1_5_",
# #         text = "two outliers removed, medNorm")
# fc_plot(handle = "log_1_5_",
#         text = "outliers removed, medNorm")
# 
# # Thesis version (fig. 1)
# fc_plot(handle = "init_",
#         text = "No outliers removed", height = 800, width = 1067,
#         top100 = FALSE)
# fc_plot(handle = "init_1_5_",
#         text = "No outliers removed, medNorm", height = 800, width = 1067,
#         top100 = FALSE)
# fc_plot(handle = "log_",
#         text = "Outliers removed", height = 800, width = 1067,
#         top100 = FALSE)


# # Checking for overlap btw. significant targets
# JYNR210602_sig <- get_signif_targets(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
#                                      group1 = 3, group2 = 1)
# 
# JYNR220601_sig <- get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
#                                      group1 = 2, group2 = 1)
# sig_1 <- rownames(JYNR210602_sig[JYNR210602_sig$sign == TRUE,])
# sig_2 <- rownames(JYNR220601_sig[JYNR220601_sig$sign == TRUE,])
# sum(!sig_1 %in% sig_2) # 20 not in SomaMouse2
# 
# sig_col <- rep("lightgrey", length(X))
# sig_col[rownames(JYNR210602_sig) %in% sig_x] <- cols[1]
# sig_col[rownames(JYNR210602_sig) %in% sig_y] <- cols[2]
# sig_col[rownames(JYNR210602_sig) %in% sig_x & rownames(JYNR210602_sig) %in% sig_y] <- cols[3]
