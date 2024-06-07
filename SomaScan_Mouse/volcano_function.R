## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel); library(scater); library(grid)

## Load data
# load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")

# Load function
source("scripts/SomaScan_Mouse/signif_targets_function.R")

# soma1 <- log_1_4_JYNR210602; soma2 <- log_1_4_JYNR220601;
# meta1 <- meta_mouse_JYNR210602; meta2 = meta_mouse_JYNR210602
# group1 <- 3; group2 <- 1; mtc <- "bh"; alpha <- 0.05; fc_cutoff <- 0.3
# log_t <- TRUE; title <- "Test"; meta_aligned = FALSE; legend.dir = "vertical"
# volcano_plot(soma1 = z_log_JYNR210602)
# volcano_plot(soma1 = z_JYNR210602)
# volcano_plot(soma1 = z_JYNR210602, log_t = FALSE)


## Wrangle and plot data
# Plotting function
volcano_plot <- function(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                         case = 3, control = 1, 
                         soma2 = soma1, meta2 = meta1,
                         group1 = case, group2 = control,
                         alpha = 0.05, mtc = "bh", fc_cutoff = 0.3,
                         log_t = TRUE, title = TRUE, legend.dir = "horizontal",
                         meta_aligned = FALSE, plt_col_scheme = "intra-study",
                         custom_col_targets = c("none", "none"), custom_col_name = "none",
                         title_vjust = -2.5){
  
  # Get log2fc, p-values, and differentially expressed targets in a df
  volc <- get_signif_targets(soma1 = soma1, meta1 = meta1, soma2 = soma2, meta2 = meta2,
                             group1 = group1, group2 = group2, alpha = alpha, mtc = mtc,
                             log_t = log_t, fc_cutoff = fc_cutoff, meta_aligned = meta_aligned)
  
  # Get number of differentially expressed targets
  n_de <- length(volc$difexpressed[!(volc$difexpressed == "NO")])
  
  # Format MTC correction for plotting
  mtc <- stringr::str_to_title(mtc)
  if (mtc == "Bh"){
    mtc <- "BH"
  }
  
  # Get adjusted alpha level
  alpha_adj <- alpha
  if (mtc == "Bonferroni"){
    alpha_adj <- alpha / dim(soma1)[2]
  }
  
  # set alpha line to the lowest accepted p-value
  if (mtc == "BH"){
    alpha_adj <- max(volc$pval[volc$pval_adj < alpha])
  }
  
  # Set plotting colours and labels
  if (group1 == 3 & group2 == 2){
    plt_col_scheme <- "DIO-NASH Vehicle"
  }
  if (group1 == 1 & group2 == 1){
    plt_col_scheme <- "Chow Vehicle"
  }
  
  if (plt_col_scheme == "intra-study"){
    if (unique(meta1$dataset) == "SomaMouse1"){
      plt_cols <- c("#53BF2E", "#FFE0C6", "#BF3D35")
    } else {
      plt_cols <- c("#53BF2E", "#C6E5FF", "#BF3D35")
    }
  }
  if (plt_col_scheme == "Chow Vehicle"){
    plt_cols <- c("#55B1FF", "#C6EAB9", "#FFA355")
  }
  if (plt_col_scheme == "DIO-NASH Vehicle"){
    plt_cols <- c("#55B1FF", "#EABEBC", "#FFA355")
  }
  
  plt_labs <- c("Downregulated", "Not significant", "Upregulated")
  
  if (length(unique(volc$difexpressed)) < 3){
    if (sum(volc$difexpressed == "UP") == 0){
      plt_cols <- plt_cols[-3]
      plt_labs <- plt_labs[-3]
    }
    if (sum(volc$difexpressed == "DOWN") == 0){
      plt_cols <- plt_cols[-1]
      plt_labs <- plt_labs[-1]
    }
  }
  
  # Add in custom colours
  if (!custom_col_targets[1] == "none"){
    custom_col_name <- paste0("_", custom_col_name)
    plt_cols <- c("blue", plt_cols)
    plt_labs <- c(custom_col_name, plt_labs)
    volc$difexpressed[volc$seqid %in% custom_col_targets] <- custom_col_name
  }
  
  # Set theme
  theme_set(theme_classic() +
              theme(plot.title = element_text(hjust = 0.5, size = 12, vjust = title_vjust),
                    legend.direction = legend.dir,
                    legend.position = "top"))
  
  # Define plot
  p <- volc %>%
    ggplot(aes(x = log2fc,
               y = -log10(pval),
               col = difexpressed)) +
    geom_hline(yintercept = -log10(alpha_adj), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = plt_cols,
                       labels = plt_labs) +
    labs(color = 'Significance', #legend_title,
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    geom_text_repel(data = subset(volc, !(difexpressed == "NO")),
                    aes(label = gene)) +
    scale_x_continuous(breaks = seq(-round(max(abs(c(min(volc$log2fc), max(volc$log2fc)))),1), round(max(abs(c(min(volc$log2fc), max(volc$log2fc)))),1), round(round(max(abs(c(min(volc$log2fc), max(volc$log2fc)))),1)/4,1))) +
    annotation_custom(grobTree(textGrob(paste0("Correction: ", mtc, ",\n",
                                               "DE targets: ", n_de, ",\n",
                                               "Alpha = ", alpha), x = 0.02, y = 0.93,
                                        hjust = 0, vjust = 1, gp = gpar(col = "black"))))
  
  if (!fc_cutoff == 0){
    p <- p +
      geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), col = "gray", linetype = 'dashed')
  }
  
  if (title == TRUE){
    p <- p + ggtitle(paste0("Group ", group1, " vs group ", group2))
  } else if (!title == FALSE){
    p <- p + ggtitle(title)
  }
  
  # Plot it
  p
}
