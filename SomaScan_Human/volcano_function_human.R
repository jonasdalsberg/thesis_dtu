## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel); library(scater); library(grid)

## Wrangle and plot data
# Plotting function
volcano_plot <- function(df1 = 1, df2 = df1, norm = "tier 1",
                         title = TRUE, legend.dir = "horizontal", title_vjust = -2.5, 
                         alpha = 0.05, fc_cutoff = 0.3, short_title = FALSE, shorter_title = TRUE,
                         frontpage = FALSE){
  
  # Extract df
  analysis <- paste0("Treatment STEP ", df1, " vs Control STEP ", df2)
  volc <- soma_human[soma_human$analysis == analysis & soma_human$normalization_method == norm,]
  
  # Get adjusted alpha level (BH MTC)
  alpha_adj <- max(volc$pvalue[volc$padj < alpha])
  
  # Set plotting colours and labels
  plt_cols <- c("#F8766D", "grey", "#00BFC4")
  col1 <- "#E0B000"
  col2 <- "#69E000"
  if (df1 == 1){
    plt_cols[2] <- "#F0D880"
  }
  if (df1 == 2){
    plt_cols[2] <- "#B4F080"
  }
  if (df1 == df2){
    title_text <- paste0("Intra-study volcano plot, SomaHuman", df1, "\nTreatment vs. placebo")
  } else {
    title_text <- paste0("Placebo-swapped volcano plot,\nSomaHuman", df1, " (treatment) vs. SomaHuman", df2,
                         " (placebo)")
  }

  
  plt_labs <- c("Downregulated", "Not significant", "Upregulated")
  
  volc$difexpressed <- "NO"
  volc$difexpressed[volc$logFC > fc_cutoff & volc$padj < alpha] <- "UP"
  volc$difexpressed[volc$logFC < -fc_cutoff & volc$padj < alpha] <- "DOWN"
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
    ggplot(aes(x = logFC,
               y = -log10(pvalue),
               col = difexpressed)) +
    geom_hline(yintercept = -log10(alpha_adj), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_colour_manual(values = plt_cols,
                       labels = plt_labs) +
    labs(colour = 'Significance', #legend_title,
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1),
                       breaks = seq(-round(max(abs(c(min(volc$logFC), max(volc$logFC)))),1),
                                    round(max(abs(c(min(volc$logFC), max(volc$logFC)))),1),
                                    round(round(max(abs(c(min(volc$logFC), max(volc$logFC)))),1)/4,1)))
  
  if (frontpage == FALSE){
    p <- p +
    annotation_custom(grobTree(textGrob(paste0("Correction: BH,\n",
                                               "DE targets: ", n_de, ",\n",
                                               "Alpha = 0.05"), x = 0.02, y = 0.93,
                                        hjust = 0, vjust = 1, gp = gpar(col = "black"))))
  }
  
  if (!fc_cutoff == 0){
    p <- p +
      geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), col = "gray", linetype = 'dashed')
  }
  
  if (shorter_title == TRUE){
    short_title <- TRUE
  }
  if (short_title == TRUE){
    if (df1 == df2){
      title_text <- paste0(str_split(title_text, ", ")[[1]][-1])
    } else {
      if (shorter_title == TRUE){
        title_text <- paste0("SomaHuman", df1, " vs. SomaHuman", df2)
      } else {
        title_text <- paste0(str_split(title_text, ",\n")[[1]][-1])
      }
    }
  }
  if (title == TRUE){
    p <- p + ggtitle(title_text)
  } else if (!title == FALSE){
    title_text <- title
    p <- p + ggtitle(title_text)
  }
  
  if (frontpage == TRUE){
    # Remove labels, numbers, titles, and change the colour scheme
    p <- p +
      labs(colour = "", x = "", y = "") +
      ggtitle("") +
      scale_colour_manual(values = c("#c51137CC", "#63154E99", "#001965CC")) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())
  }
  
  # Plot it
  p
}
