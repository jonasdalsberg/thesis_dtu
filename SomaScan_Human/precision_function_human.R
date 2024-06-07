## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel); library(scater); library(grid)

## Functions
get_signif_targets <- function(df1 = 1, df2 = df1, norm = "tier 1",
                               fc_cutoff = 0.3, alpha = 0.05, top100 = FALSE){

  # Extract df
  analysis <- paste0("Treatment STEP ", df1, " vs Control STEP ", df2)
  volc <- soma_human[soma_human$analysis == analysis & soma_human$normalization_method == norm,]
  
  if (top100 == TRUE){
    volc <- volc[order(abs(volc$logFC), decreasing = TRUE),]
    volc <- volc[abs(volc$logFC) > fc_cutoff & volc$padj < alpha,]
    volc <- volc[1:100, colnames(volc) == "variable"]
    return(volc)
  }
  
  # Get list of significant (DE) target
  volc$difexpressed <- "NO"
  volc$difexpressed[volc$logFC > fc_cutoff & volc$padj < alpha] <- "UP"
  volc$difexpressed[volc$logFC < -fc_cutoff & volc$padj < alpha] <- "DOWN"
  volc$sign <- FALSE
  volc$sign[!volc$difexpressed == "NO"] <- TRUE
  sig_list <- volc$variable[volc$sign == TRUE]
  
  return(sig_list)
}

## Wrangle and plot data
# Plotting function
compute_precision <- function(df1 = 1, norm = "tier 1",
                         alpha = 0.05, fc_cutoff = 0.3, extra_info = FALSE){
  
  if (df1 == 1){
    df2 <- 2
  } else {
    df2 <- 1
  }
  
  intra_list_1 <- get_signif_targets(df1 = df1, df2 = df1, fc_cutoff = fc_cutoff,
                                     alpha = alpha, norm = norm)
  intra_list_2 <- get_signif_targets(df1 = df2, df2 = df2, fc_cutoff = fc_cutoff,
                                     alpha = alpha, norm = norm)
  swap_list <- get_signif_targets(df1 = df1, df2 = df2, fc_cutoff = fc_cutoff,
                                  alpha = alpha, norm = norm)
  top100_1 <- get_signif_targets(df1 = df1, df2 = df1, fc_cutoff = fc_cutoff,
                                 alpha = alpha, top100 = TRUE, norm = norm)
  top100_2 <- get_signif_targets(df1 = df2, df2 = df2, fc_cutoff = fc_cutoff,
                                 alpha = alpha, top100 = TRUE, norm = norm)
  
  intra_no <- length(intra_list_1)
  swap_no <- length(swap_list)
  precision <- sum(intra_list_1 %in% swap_list) / length(swap_list)
  overlap <- sum(intra_list_1 %in% intra_list_2) / length(intra_list_1)
  top100 <- sum(top100_1 %in% top100_2)
  
  if (norm == "raw"){
    norm_txt <- "Unnorm"
  }
  if (norm == "tier 1"){
    norm_txt <- "Base"
  }
  if (norm == "ANML"){
    norm_txt <- norm
  }

  dge_results <- data.frame(`intra-study` = intra_no,
                            `placebo-swapped` = swap_no,
                            precision = precision,
                            overlap = overlap,
                            top100 = top100,
                            dataset = paste0("SomaHuman", df1),
                            norm = norm_txt)
  if (extra_info == TRUE){
    dge_results$fc_cutoff <- fc_cutoff
    dge_results$setting <- paste0(df1, " vs. ", df2)
    dge_results$analysis <- paste0(norm_txt, " ", df1)
    dge_results$analysis_exp <- paste0(norm_txt, "_", df1, "_", fc_cutoff)
  }
  rownames(dge_results) <- paste0(norm_txt, "_", df1, "_", fc_cutoff)
  return(dge_results)
}


# Auxilliary function
merge_dfs <- function(df1, df2){
  ex <- merge(t(df1), t(df2), by = "row.names", all = TRUE)
  rownames(ex) <- ex[,1]
  ex <- ex[-1]
  ex <- as.data.frame(t(ex))
  return (ex)
}
