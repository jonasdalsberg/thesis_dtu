## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel); library(scater); library(grid)

## Wrangle and plot data
# Plotting function
volcano_plot <- function(df1 = 1, df2 = df1, type = "intra",
                         title = TRUE, legend.dir = "horizontal", 
                         plt_col_scheme = "intra-study", title_vjust = -2.5, 
                         alpha = 0.05, fc_cutoff = 0.3, short_title = FALSE,
                         lab = TRUE){
  
  # Extract df
  df_name <- paste0("dge_t", df1, "_", type)
  if (!df2 == df1){
    df_name <- paste0(df_name, "_", df2)
  }
  dge_df <- get(df_name)
  
  # Remove padj NAs
  volc <- dge_df[!is.na(dge_df$padj),]
  
  # Get adjusted alpha level (BH MTC)
  alpha_adj <- max(volc$pvalue[volc$padj < alpha])
  
  # Set plotting colours and labels
  type_split <- str_split(type, "_")[[1]]
  plt_cols <- c("#53BF2E", "grey", "#BF3D35")
  if (type == "intra" || type == "swap"){
    if (df1 == 1){
      plt_cols[2] <- "#BBB9EA"
    }
    if (df1 == 2){
      plt_cols[2] <- "#D4B9EA"
    }
    if (df1 == 3){
      plt_cols[2] <- "#EAB9E3"
    }
    if (type == "intra"){
      title_text <- paste0("Intra-study volcano plot, TransMouse", df1, "\nDIO-NASH vehicle vs. chow vehicle")
    } else {
      title_text <- paste0("Placebo-swapped volcano plot, TransMouse", df1, " vs. TransMouse", df2,
                           "\nDIO-NASH vehicle vs. chow vehicle")
    }
  } else {
    if (type_split[1] == "inter"){
      
      if (df1 == 1){
        plt_cols[3] <- "#332EBF"
      }
      if (df1 == 2){
        plt_cols[3] <- "#7D2EBF"
      }
      if (df1 == 3){
        plt_cols[3] <- "#BF2EAC"
      }
      if (df2 == 1){
        plt_cols[1] <- "#332EBF"
      }
      if (df2 == 2){
        plt_cols[1] <- "#7D2EBF"
      }
      if (df2 == 3){
        plt_cols[1] <- "#BF2EAC"
      }
      if (type_split[2] == "chow"){
        plt_cols[2] <- "#C6EAB9"
        title_text <- paste0("Inter-study volcano plot, TransMouse", df1,
                             " vs. TransMouse", df2, "\nChow vehicle vs. chow vehicle")
      } else {
        plt_cols[2] <- "#EABEBC"
        title_text <- paste0("Inter-study volcano plot, TransMouse", df1,
                             " vs. TransMouse", df2, "\nDIO-NASH vehicle vs DIO-NASH vehicle")
      }
    }
  }

  
  plt_labs <- c("Downregulated", "Not significant", "Upregulated")
  
  volc$difexpressed <- "NO"
  volc$difexpressed[volc$log2FoldChange > fc_cutoff & volc$padj < alpha] <- "UP"
  volc$difexpressed[volc$log2FoldChange < -fc_cutoff & volc$padj < alpha] <- "DOWN"
  volc$sign <- FALSE
  volc$sign[!volc$difexpressed == "NO"] <- TRUE
  
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
  
  # Get number of differentially expressed targets
  n_de <- sum(!volc$difexpressed == "NO")
  
  # Set theme
  theme_set(theme_classic() +
              theme(plot.title = element_text(hjust = 0.5, size = 12, vjust = title_vjust),
                    legend.direction = legend.dir,
                    legend.position = "top"))
  
  # Define plot
  p <- volc %>%
    ggplot(aes(x = log2FoldChange,
               y = -log10(pvalue),
               col = difexpressed)) +
    geom_hline(yintercept = -log10(alpha_adj), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = plt_cols,
                       labels = plt_labs) +
    labs(color = 'Significance', #legend_title,
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    scale_x_continuous(breaks = seq(-round(max(abs(c(min(volc$log2FoldChange), max(volc$log2FoldChange)))),1),
                                    round(max(abs(c(min(volc$log2FoldChange), max(volc$log2FoldChange)))),1),
                                    round(round(max(abs(c(min(volc$log2FoldChange), max(volc$log2FoldChange)))),1)/4,1))) +
    annotation_custom(grobTree(textGrob(paste0("Correction: BH,\n",
                                               "DE targets: ", n_de, ",\n",
                                               "Alpha = 0.05"), x = 0.02, y = 0.93,
                                        hjust = 0, vjust = 1, gp = gpar(col = "black"))))
  
  if (lab == TRUE){
    p <- p + geom_text_repel(data = subset(volc, !(difexpressed == "NO")),
                                 aes(label = name))
  }
  
  if (!fc_cutoff == 0){
    p <- p +
      geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), col = "gray", linetype = 'dashed')
  }
  
  if (short_title == TRUE){
    title_text <- paste0(str_split(title_text, ", ")[[1]][-1], collapse = ", ")
  }
  if (title == TRUE){
    p <- p + ggtitle(title_text)
  } else if (!title == FALSE){
    title_text <- title
    p <- p + ggtitle(title_text)
  }
  
  # Plot it
  p
}
# 
# volcano_plot()
# volcano_plot(df2 = 2, type = "inter_nash")
# volcano_plot(df2 = 3, type = "inter_chow")
# volcano_plot(df1 = 3, df2 = 2, type = "swap")
#
# 
# test_string <- c("Text goes, here with, my girl")
# paste0(str_split(test_string, ", ")[[1]][-1], collapse = ", ")
# 
# text <- str_split(test_string, ",")[[1]][-1]
# text
# paste(text, collapse = ",")
