## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel); library(scater); library(grid)


## Functions
get_signif_targets <- function(df1 = 1, df2 = df1, type = "intra",
                               fc_cutoff = 0.3, alpha = 0.05, top100 = FALSE){

  # Extract df
  df_name <- paste0("dge_t", df1, "_", type)
  if (!df2 == df1){
    df_name <- paste0(df_name, "_", df2)
  }
  dge_df <- get(df_name)
  
  # Remove padj NAs
  volc <- dge_df[!is.na(dge_df$padj),]
  
  if (top100 == TRUE){
    volc <- volc[order(abs(volc$log2FoldChange), decreasing = TRUE),]
    volc <- volc[abs(volc$log2FoldChange) > fc_cutoff & volc$padj < alpha,]
    volc <- volc[1:100, colnames(volc) == "identifier"]
    return(volc)
  }
  
  # Get list of significant (DE) target
  volc$difexpressed <- "NO"
  volc$difexpressed[volc$log2FoldChange > fc_cutoff & volc$padj < alpha] <- "UP"
  volc$difexpressed[volc$log2FoldChange < -fc_cutoff & volc$padj < alpha] <- "DOWN"
  volc$sign <- FALSE
  volc$sign[!volc$difexpressed == "NO"] <- TRUE
  sig_list <- volc$identifier[volc$sign == TRUE]
  
  return(sig_list)
}

## Wrangle and plot data
# Plotting function
compute_precision <- function(df1 = 1, df2 = 2,
                         alpha = 0.05, fc_cutoff = 0.3){
  intra_list_1 <- get_signif_targets(df1 = df1, df2 = df1, type = "intra", fc_cutoff = fc_cutoff, alpha = alpha)
  intra_list_2 <- get_signif_targets(df1 = df2, df2 = df2, type = "intra", fc_cutoff = fc_cutoff, alpha = alpha)
  swap_list <- get_signif_targets(df1 = df1, df2 = df2, type = "swap", fc_cutoff = fc_cutoff, alpha = alpha)
  top100_1 <- get_signif_targets(df1 = df1, df2 = df1, type = "intra",
                                    fc_cutoff = fc_cutoff, alpha = alpha, top100 = TRUE)
  top100_2 <- get_signif_targets(df1 = df2, df2 = df2, type = "intra",
                                      fc_cutoff = fc_cutoff, alpha = alpha, top100 = TRUE)
  
  intra_no <- length(intra_list_1)
  swap_no <- length(swap_list)
  precision <- sum(intra_list_1 %in% swap_list) / length(swap_list)
  overlap <- sum(intra_list_1 %in% intra_list_2) / length(intra_list_1)
  top100 <- sum(top100_1 %in% top100_2)

  dge_results <- data.frame(`intra-study` = intra_no,
                            `placebo-swapped` = swap_no,
                            precision = precision,
                            overlap = overlap,
                            top100 = top100,
                            dataset = paste0("TransMouse", df1),
                            setting = paste0(df1, " vs. ", df2))
  rownames(dge_results) <- paste0("TransMouse", df1, "_vs_TransMouse", df2)
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
