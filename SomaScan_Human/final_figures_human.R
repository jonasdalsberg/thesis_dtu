### Initialize document
# Load libraries
library(tidyverse)

# Set working directory
setwd("~/JDJG/clean/")

# Load functions
source("scripts/SomaScan_Human/volcano_function_human.R")
source("scripts/SomaScan_Human/precision_function_human.R")

# Load data
load("~/NNEDL/masterprojectxjdjg/curated/Human_studies/soma_human.RData")


################################################################################
################################################################################
################################# MAIN FIGURES #################################
################################################################################
################################################################################


################################################################################
#################### FIGURE 7: Volcano, precision, overlap #####################
################################################################################
for (generate_figure in "FIGURE_7"){

  
  ### FIGURE 7A: Intra-study volcano plots
  myplots <- list(volcano_plot(df1 = 1, short_title = TRUE),
                  volcano_plot(df1 = 2, short_title = TRUE))
  g <- ggarrange(plotlist = myplots,
                 nrow = 1, ncol = 2,
                 hjust = -2)
  fig7a <- annotate_figure(g, top = text_grob("Intra-study volcano plots",
                                              vjust = 1.5, face = "bold"))
  
  
  ### FIGURE 7B: Placebo-swapped volcano plots
  myplots <- list(volcano_plot(df1 = 1, df2 = 2, short_title = TRUE),
                  volcano_plot(df1 = 2, df2 = 1, short_title = TRUE))
  g <- ggarrange(plotlist = myplots,
                 nrow = 1, ncol = 2,
                 hjust = -2)
  fig7b <- annotate_figure(g, top = text_grob("Placebo-swapped volcano plots",
                                              vjust = 1.5, face = "bold"))
  
  
  ### FIGURE 7C+7D
  ## Common setup
  # Get precision and overlap data
  cnt <- 0
  for (norm in c("raw", "tier 1", "ANML")){
    for (df in c(1,2)){
      for (fc_cutoff in c(0.3)){
        cnt <- cnt + 1
        pre <- compute_precision(df1 = df, fc_cutoff = fc_cutoff, norm = norm)
        if (cnt == 1){
          precision_data <- pre
        } else {
          precision_data <- merge_dfs(precision_data, pre)
        }
      }
    }
  }
  
  for (i in c(2,4:7)){precision_data[,i] <- as.numeric(precision_data[,i])}
  precision_data$norm <- factor(precision_data$norm,
                                levels = c("Unnorm", "Base",
                                           "ANML"))
  
  # Define colours
  col1 <- "#E0B000"
  col2 <- "#69E000"
  
  # Compute positions for dashed lines
  x_start <- seq(0.5,3,0.5); x_end <- x_start + 0.5
  
  
  ## FIGURE 7C
  # Compute positions for dashed lines
  y_start_1 <- c()
  i_ <- 0.5
  for (i in c(1,3,5)){
    pre1 <- precision_data$precision[i]
    pre2 <- precision_data$precision[i+1]
    med_pre <- median(c(pre1,pre2))
    y_start_1 <- c(y_start_1, rep(med_pre,2))
  }
  y_end_1 <- y_start_1
  
  # Plot
  fig7c_1 <- precision_data %>% 
    ggplot(mapping = aes(x = norm,
                         y = precision,
                         fill = dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1, vjust = 1),
          legend.direction = "horizontal") +
    ylim(c(0,1)) +
    scale_fill_manual(values = rep(c(col1, col2),3)) +
    xlab("Normalisation method") + ylab("Precision") + labs(fill = "Dataset") +
    geom_segment(aes(x = x_start, xend = x_end,
                     y = y_start_1,
                     yend = y_end_1),
                 linetype = "dashed", color = rep("#A5C800", 6))
  
  # Compute positions for dashed lines
  y_end_2 <- y_start_2 <- precision_data$intra.study
  
  # Plot DE targets (placebo-swapped with intra-study as a reference)
  fig7c_2 <- precision_data %>%
    ggplot(mapping = aes(x = norm,
                         y = placebo.swapped,
                         fill = dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1, vjust = 1),
          legend.direction = "vertical") +
    scale_fill_manual(values = rep(c(col1, col2),3)) +
    xlab("Normalisation method") + ylab("DE targets") + labs(fill = "Dataset") +
    geom_segment(aes(x = x_start, xend = x_end,
                     y = y_start_2,
                     yend = y_end_2),
                 linetype = "dashed", color = rep(c(col1, col2),3))
  
  legend <- get_legend(fig7c_1)
  g <- ggarrange(fig7c_1,fig7c_2,
                 ncol = 2,
                 common.legend = TRUE,
                 legend.grob = legend,
                 hjust = -2)
  title <- "Evaluation of swapping feasibility"
  fig7c <- annotate_figure(g, top = text_grob(title, vjust = 1, face = "bold"))
  
  
  ## FIGURE 7D
  # Compute positions for dashed lines
  y_start_1 <- c(); i_ <- 0.5
  for (i in c(1,3,5)){
    overlap1 <- precision_data$overlap[i]
    overlap2 <- precision_data$overlap[i+1]
    med_overlap <- median(c(overlap1,overlap2))
    y_start_1 <- c(y_start_1, rep(med_overlap,2))
  }
  y_end_1 <- y_start_1
  
  # Plot Fig. 7d_2
  fig7d <- precision_data %>% 
    ggplot(mapping = aes(x = norm,
                         y = overlap,
                         fill = dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1, vjust = 1),
          legend.direction = "horizontal") +
    ylim(c(0,1)) +
    scale_fill_manual(values = rep(c(col1, col2),3)) +
    xlab("Normalisation method") + ylab("Total overlap") + labs(fill = "Dataset") +
    geom_segment(aes(x = x_start, xend = x_end,
                     y = y_start_1,
                     yend = y_end_1),
                 linetype = "dashed", color = rep("#A5C800", 6))
  
  legend <- get_legend(fig7d)
  g <- ggarrange(fig7d,
                 ncol = 1,
                 common.legend = TRUE,
                 legend.grob = legend,
                 hjust = -2)
  title <- "Inter-study overlap"
  fig7d <- annotate_figure(g, top = text_grob(title, vjust = 1, face = "bold"))
  
  # will probably change this figure to be part mouse transcriptomics and part human proteomics
  # as I do not feel that the conclusions are strong enough for a full figure
  
  ## Put it all together
  fig7 <- ggarrange(fig7a,fig7b,fig7c,fig7d,
                    nrow = 4, heights = c(1.5,1.5,1,1),
                    labels = c("A", "B", "C", "D"))
  
  ## Save it
  fact <- 1.5
  ggsave("figures/fig7.png", plot = fig7, width = fact*(210-20), height = fact*(297-30),
         units = "mm", dpi = 185, bg = "white", limitsize = FALSE)
}



################################################################################
############################ Supplementary figures #############################
################################################################################

################################################################################
## Supplementary figure 19: Human data (intra-study volcano)
for (generate_figure in "sup_fig_19"){
  
  # Unnormalised data
  myplots <- list(volcano_plot(df1 = 1, short_title = TRUE, norm = "raw"),
                  volcano_plot(df1 = 2, short_title = TRUE, norm = "raw"))
  g <- ggarrange(plotlist = myplots,
                 nrow = 1, ncol = 2,
                 hjust = -2)
  sup_fig_19a <- annotate_figure(g, top = text_grob("Intra-study volcano plots\nUnnormalised data",
                                              vjust = 1.5, face = "bold"))
  
  # Baseline data
  myplots <- list(volcano_plot(df1 = 1, short_title = TRUE, norm = "tier 1"),
                  volcano_plot(df1 = 2, short_title = TRUE, norm = "tier 1"))
  g <- ggarrange(plotlist = myplots,
                 nrow = 1, ncol = 2,
                 hjust = -2)
  sup_fig_19b <- annotate_figure(g, top = text_grob("Intra-study volcano plots\nBaseline data",
                                                    vjust = 1.5, face = "bold"))
  
  # ANML data
  myplots <- list(volcano_plot(df1 = 1, short_title = TRUE, norm = "ANML"),
                  volcano_plot(df1 = 2, short_title = TRUE, norm = "ANML"))
  g <- ggarrange(plotlist = myplots,
                 nrow = 1, ncol = 2,
                 hjust = -2)
  sup_fig_19c <- annotate_figure(g, top = text_grob("Intra-study volcano plots\nANML data",
                                                    vjust = 1.5, face = "bold"))
  
  sup_fig_19 <- ggarrange(sup_fig_19a, sup_fig_19b, sup_fig_19c,
                          nrow = 3, ncol = 1,
                          labels = c("A", "B", "C"))
  
  fact <- 1.5
  ggsave("figures/supfig19.png", plot = sup_fig_19, width = fact*(210-20), height = fact*(297-30-35),
         units = "mm", dpi = 185, bg = "white")
}

################################################################################
## Supplementary figure 20: Human data (inter-study, chow)
for (generate_figure in "sup_fig_20"){
  
  # Unnormalised data
  myplots <- list(volcano_plot(df1 = 1, df2 = 2, shorter_title = TRUE, norm = "raw"),
                  volcano_plot(df1 = 2, df2 = 1, shorter_title = TRUE, norm = "raw"))
  g <- ggarrange(plotlist = myplots,
                 nrow = 1, ncol = 2,
                 hjust = -2)
  sup_fig_20a <- annotate_figure(g, top = text_grob("Placebo-swapped volcano plots\nUnnormalised data",
                                                    vjust = 1.5, face = "bold"))
  
  # Baseline data
  myplots <- list(volcano_plot(df1 = 1, df2 = 2, shorter_title = TRUE, norm = "tier 1"),
                  volcano_plot(df1 = 2, df2 = 1, shorter_title = TRUE, norm = "tier 1"))
  g <- ggarrange(plotlist = myplots,
                 nrow = 1, ncol = 2,
                 hjust = -2)
  sup_fig_20b <- annotate_figure(g, top = text_grob("Placebo-swapped volcano plots\nBaseline data",
                                                    vjust = 1.5, face = "bold"))
  
  # ANML data
  myplots <- list(volcano_plot(df1 = 1, df2 = 2, shorter_title = TRUE, norm = "ANML"),
                  volcano_plot(df1 = 2, df2 = 1, shorter_title = TRUE, norm = "ANML"))
  g <- ggarrange(plotlist = myplots,
                 nrow = 1, ncol = 2,
                 hjust = -2)
  sup_fig_20c <- annotate_figure(g, top = text_grob("Placebo-swapped volcano plots\nANML data",
                                                    vjust = 1.5, face = "bold"))
  
  sup_fig_20 <- ggarrange(sup_fig_20a, sup_fig_20b, sup_fig_20c,
                          nrow = 3, ncol = 1,
                          labels = c("A", "B", "C"))
  
  fact <- 1.5
  ggsave("figures/supfig20.png", plot = sup_fig_20, width = fact*(210-20), height = fact*(297-30-35),
         units = "mm", dpi = 185, bg = "white")
}

################################################################################
## Supplementary figure 21: Cross-normalisation overlap
for (generate_figure in "sup_fig_21"){
  
  ## Pre-processing the data
  for (preprocess in "now"){
    base_1 <- get_signif_targets(df1 = 1)
    base_2 <- get_signif_targets(df1 = 2)
    anml_1 <- get_signif_targets(df1 = 1, norm = "ANML")
    anml_2 <- get_signif_targets(df1 = 2, norm = "ANML")
    base_swap_1 <- get_signif_targets(df1 = 1, df2 = 2)
    base_swap_2 <- get_signif_targets(df1 = 2, df2 = 1)
    anml_swap_1 <- get_signif_targets(df1 = 1, df2 = 2, norm = "ANML")
    anml_swap_2 <- get_signif_targets(df1 = 2, df2 = 1, norm = "ANML")
    
    base_1_no <- length(base_1)
    base_2_no <- length(base_2)
    anml_1_no <- length(anml_1)
    anml_2_no <- length(anml_2)
    
    base_swap_1_no <- length(base_swap_1)
    base_swap_2_no <- length(base_swap_2)
    anml_swap_1_no <- length(anml_swap_1)
    anml_swap_2_no <- length(anml_swap_2)
    
    base_1_overlap <- sum(base_1 %in% anml_1) / length(base_1)
    base_2_overlap <- sum(base_2 %in% anml_2) / length(base_2)
    anml_1_overlap <- sum(base_1 %in% anml_1) / length(anml_1)
    anml_2_overlap <- sum(base_2 %in% anml_2) / length(anml_2)
    
    base_swap_1_overlap <- sum(base_swap_1 %in% anml_swap_1) / length(base_swap_1)
    base_swap_2_overlap <- sum(base_swap_2 %in% anml_swap_2) / length(base_swap_2)
    anml_swap_1_overlap <- sum(base_swap_1 %in% anml_swap_1) / length(anml_swap_1)
    anml_swap_2_overlap <- sum(base_swap_2 %in% anml_swap_2) / length(anml_swap_2)
    
    norm_overlap <- data.frame(dataset = rep(c("SomaHuman1", "SomaHuman1", "SomaHuman2", "SomaHuman2"), 2),
                               norm = factor(rep(c("Base", "ANML"), 4),
                                             levels = c("Base", "ANML")),
                               sig.tar = c(base_1_no, anml_1_no,
                                           base_2_no, anml_2_no,
                                           base_swap_1_no, anml_swap_1_no,
                                           base_swap_2_no, anml_swap_2_no),
                               overlap = c(base_1_overlap, anml_1_overlap,
                                           base_2_overlap, anml_2_overlap,
                                           base_swap_1_overlap, anml_swap_1_overlap,
                                           base_swap_2_overlap, anml_swap_2_overlap),
                               analysis = c(rep("Intra-study", 4),
                                            rep("Placebo-swapped", 4)),
                               setting = c(rep("Intra-study 1",2),
                                           rep("Intra-study 2", 2),
                                           rep("Placebo-swapped 1", 2),
                                           rep("Placebo-swapped 2", 2)))
  }
  
  
  ## Plotting the data
  # Compute positions for dashed lines
  x_start <- c(0,0,1,1,2,2,3,3) + 0.5; x_end <- c(1,1,2,2,3,3,4,4) + 0.5
  y_start_1 <- c(); i_ <- 0.5
  for (i in c(1,3,5,7)){
    med_overlap <- median(norm_overlap$overlap[i:(i+1)])
    y_start_1 <- c(y_start_1, rep(med_overlap,2))
  }
  y_end_1 <- y_start_1
  
  # Plot
  sup_fig_21a <- norm_overlap %>% 
    ggplot(mapping = aes(x = setting,
                         y = overlap,
                         fill = dataset,
                         col = norm)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 30,
                                     hjust = 1, vjust = 1),
          legend.direction = "horizontal") +
    ylim(c(0,1)) +
    scale_fill_manual(values = rep(c("#E0B000", "#69E000"),3)) +
    scale_colour_manual(values = rep(c("#0030E0", "#7700E0"),3)) +
    xlab("Analysis") + ylab("Overlap") + labs(fill = "Dataset ", col = "Normalisation method     ") +
    geom_segment(aes(x = x_start, xend = x_end,
                     y = y_start_1,
                     yend = y_end_1),
                 linetype = "dashed", color = rep("#A5C800", 8))
  
  
  
  # Compute positions for dashed lines
  y_start_2 <- c(); i_ <- 0.5
  for (i in c(1,3,5,7)){
    med_de <- median(norm_overlap$sig.tar[i:(i+1)])
    y_start_2 <- c(y_start_2, rep(med_de,2))
  }
  y_end_2 <- y_start_2
  
  # Plot DE targets (placebo-swapped with intra-study as a reference)
  sup_fig_21b <- norm_overlap %>% 
    ggplot(mapping = aes(x = setting,
                         y = sig.tar,
                         fill = dataset,
                         col = norm)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 30,
                                     hjust = 1, vjust = 1),
          legend.direction = "horizontal") +
    scale_fill_manual(values = rep(c("#E0B000", "#69E000"),3)) +
    scale_colour_manual(values = rep(c("#0030E0", "#7700E0"),3)) +
    xlab("Analysis") + ylab("DE targets") + labs(fill = "Dataset", col = "Normalisation method") +
    geom_segment(aes(x = x_start, xend = x_end,
                     y = y_start_2,
                     yend = y_end_2),
                 linetype = "dashed", color = rep("#A5C800", 8))
  
  g <- ggarrange(sup_fig_21a,sup_fig_21b,
                 ncol = 2,
                 common.legend = TRUE,
                 hjust = -2)
  title <- "Cross-normalisation overlap"
  sup_fig_21 <- annotate_figure(g, top = text_grob(title, vjust = 1, face = "bold"))
  fact <- 1.5
  ggsave("figures/supfig21.png", plot = sup_fig_21, width = fact*(210-20), height = fact*(297-30-220),
         units = "mm", dpi = 185, bg = "white")
  
  # The "Normalisation method" legend label is cut off by its own content.
  # Fix it.
  # A quickfix would be to shorten the legend label. Maybe it might even work to switch the position of the legends.
  
  # cowplot::save_plot("figures/supfig21.png", plot = g, base_width = fact*(210-20) / 25.4, base_height = fact*(297-30-150) / 25.4)
  # ?cowplot::save_plot()
}

################################################################################
## Supplementary figure 22: Optimising FC for different normalisation methods
for (generate_figure in "sup_fig_22"){
  # Generate inter-study volcano plots (lean chow)
  # Get precision and overlap data
  cnt <- 0
  for (norm in c("raw", "tier 1", "ANML")){
    for (df1 in c(1,2)){
      for (fc_cutoff in c(0,0.1,0.2,0.3,0.4,0.5)){
        cnt <- cnt + 1
        pre <- compute_precision(df1 = df1, fc_cutoff = fc_cutoff, norm = norm,
                                 extra_info = TRUE)
        if (cnt == 1){
          precision_data <- pre
        } else {
          precision_data <- merge_dfs(precision_data, pre)
        }
      }
    }
  }
  
  for (i in c(5,7,8,9,11)){precision_data[,i] <- as.numeric(precision_data[,i])}
  precision_data$analysis <- factor(precision_data$analysis,
                                    levels = c("Unnorm 1", "Unnorm 2",
                                               "Base 1", "Base 2",
                                               "ANML 1", "ANML 2"))
  # Remove unnormalised data
  precision_data <- precision_data[!precision_data$analysis %in% c("Unnorm 1", "Unnorm 2"),]
  
  # Create extra variable
  precision_data$analysis_fc <- factor(paste0(precision_data$norm, " ", precision_data$fc_cutoff),
                                       levels = unique(paste0(precision_data$norm, " ", precision_data$fc_cutoff)))
  
  
  
  ## Supplementary figure 22a
  # Compute positions for dashed lines
  x_start <- unlist(lapply(1:12, function(x) rep(x,2)))-0.5; x_end <- x_start + 1
  base_meds <- anml_meds <- c()
  for (i in 1:6){
    j <- i + 6; k <- j + 6; l <- k + 6;
    med_pre_base <- median(precision_data$precision[c(i,j)])
    med_pre_anml <- median(precision_data$precision[c(k,l)])
    base_meds <- c(base_meds, med_pre_base, med_pre_base)
    anml_meds <- c(anml_meds, med_pre_anml, med_pre_anml)
  }
  y_start_1 <- y_end_1 <- c(base_meds, anml_meds)
  
  # Plot
  sup_fig_22a_1 <- precision_data %>% 
    ggplot(mapping = aes(x = analysis_fc,
                         y = precision,
                         fill = dataset,
                         col = norm)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 30,
                                     hjust = 1, vjust = 1),
          legend.direction = "horizontal") +
    ylim(c(0,1)) +
    scale_fill_manual(values = rep(c("#E0B000", "#69E000"),3)) +
    xlab("Analysis") + ylab("Precision") + labs(fill = "Dataset  ", col = "Normalisation method     ") +
    geom_segment(aes(x = x_start, xend = x_end,
                     y = y_start_1,
                     yend = y_end_1),
                 linetype = "dashed", color = rep("#A5C800", 24))
  
  # Compute positions for dashed lines
  x_start <- unlist(lapply(1:12, function(x) rep(x,2)))-0.5; x_end <- x_start + 1
  base_meds <- anml_meds <- c()
  for (i in 1:6){
    j <- i + 6; k <- j + 6; l <- k + 6;
    med_ov_base <- median(precision_data$overlap[c(i,j)])
    med_ov_anml <- median(precision_data$overlap[c(k,l)])
    base_meds <- c(base_meds, med_ov_base, med_ov_base)
    anml_meds <- c(anml_meds, med_ov_anml, med_ov_anml)
  }
  y_start_2 <- y_end_2 <- c(base_meds, anml_meds)
  
  # Plot DE targets (placebo-swapped with intra-study as a reference)
  sup_fig_22a_2 <- precision_data %>%
    ggplot(mapping = aes(x = analysis_fc,
                         y = overlap,
                         fill = dataset,
                         col = norm)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 30,
                                     hjust = 1, vjust = 1),
          legend.direction = "horizontal") +
    scale_fill_manual(values = rep(c("#E0B000", "#69E000"),3)) +
    xlab("Analysis") + ylab("Inter-study overlap") + labs(fill = "Dataset  ", col = "Normalisation method     ") +
    geom_segment(aes(x = x_start, xend = x_end,
                     y = y_start_2,
                     yend = y_end_2),
                 linetype = "dashed", color = rep("#A5C800", 24))
  
  legend <- get_legend(sup_fig_22a_2)
  g <- ggarrange(sup_fig_22a_1,sup_fig_22a_2,
                 ncol = 2,
                 common.legend = TRUE,
   #              legend.grob = legend,
                 hjust = -2)
  title <- "FC thresholding for baseline vs. ANML, inter-study comparison"
  sup_fig_22a <- annotate_figure(g, top = text_grob(title, vjust = 0.5, face = "bold"))
  
  
  sup_fig_22b_1 <- precision_data %>% 
    ggplot(mapping = aes(x = analysis,
                         y = intra.study,
                         fill = dataset,
                         col = fc_cutoff)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 30,
                                     hjust = 1, vjust = 1),
          legend.direction = "horizontal") +
    scale_fill_manual(values = rep(c("#E0B000", "#69E000"),3)) +
    xlab("Analysis") + ylab("DE targets (intra-study)") + labs(fill = "Dataset  ", col = "log2(FC) threshold    ")
  
  # Plot DE targets (placebo-swapped with intra-study as a reference)
  sup_fig_22b_2 <- precision_data %>%
    ggplot(mapping = aes(x = analysis,
                         y = placebo.swapped,
                         fill = dataset,
                         col = fc_cutoff)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 30,
                                     hjust = 1, vjust = 1),
          legend.direction = "vertical") +
    scale_fill_manual(values = rep(c("#E0B000", "#69E000"),3)) +
    xlab("Analysis") + ylab("DE targets (post-swap)") + labs(fill = "Dataset  ", col = "log2(FC) threshold    ")
  
  g <- ggarrange(sup_fig_22b_1,sup_fig_22b_2,
                 ncol = 2,
                 common.legend = TRUE,
                 hjust = -2)
  title <- "FC thresholding for baseline vs. ANML, DE targets"
  sup_fig_22b <- annotate_figure(g, top = text_grob(title, vjust = 0.5, face = "bold"))  
  
  sup_fig_22 <- ggarrange(plotlist = list(sup_fig_22a, sup_fig_22b),
                          legend = "none",
                          nrow = 2, ncol = 1,
                          labels = c("A", "B"))
  
  fact <- 1.5
  ggsave("figures/supfig22.png", plot = sup_fig_22, width = fact*(210-20), height = fact*(297-30-120),
         units = "mm", dpi = 185, bg = "white")
}

