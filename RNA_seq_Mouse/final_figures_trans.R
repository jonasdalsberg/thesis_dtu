### Initialize document
# Load libraries
library(tidyverse); library(DESeq2)

# Set working directory
setwd("~/JDJG/clean/")

# Load functions
source("scripts/RNA_seq_Mouse/volcano_function_trans.R")
source("scripts/RNA_seq_Mouse/precision_function_trans.R")

# Load data
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/trans_dge.RData")


################################################################################
################################################################################
################################# MAIN FIGURES #################################
################################################################################
################################################################################


################################################################################
############# FIGURE 6: Volcano, FC:FC (?), precision, overlap (?) ############# 
################################################################################
for (generate_figure in "FIGURE_6"){
  
  # Generate intra-study volcano plots
  myplots <- list(volcano_plot(short_title = TRUE, lab = FALSE),
                  volcano_plot(df1 = 2, short_title = TRUE, lab = FALSE),
                  volcano_plot(df1 = 3, short_title = TRUE, lab = FALSE))
  g <- ggarrange(plotlist = myplots,
                 nrow = 1, ncol = 3,
                 hjust = -2,
                 common.legend = TRUE)
  fig6a <- annotate_figure(g, top = text_grob("Intra-study volcano plots",
                                                   vjust = 1.5, face = "bold"))
  
  # Generate inter-study volcano plots (lean chow)
  myplots <- list(volcano_plot(df2 = 2, type = "inter_chow", short_title = TRUE, lab = FALSE),
                  volcano_plot(df2 = 3, type = "inter_chow", short_title = TRUE, lab = FALSE),
                  volcano_plot(df1 = 2, df2 = 3, type = "inter_chow", short_title = TRUE, lab = FALSE))
  g <- ggarrange(plotlist = myplots,
                 nrow = 1, ncol = 3,
                 hjust = -2,
                 common.legend = TRUE)
  fig6b <- annotate_figure(g, top = text_grob("Inter-study volcano plots",
                                                   vjust = 1.5, face = "bold"))
  
  # Get precision and overlap data
  d <- list()
  for (df1 in c(1,2,3)){
    for (df2 in c(1,2,3)){
      if (df1 == df2){
        next
      }
      d <- append(d, list(compute_precision(df1 = df1, df2 = df2, fc_cutoff = 0.3)))
    }
  }
  precision_data <- merge_dfs(merge_dfs(merge_dfs(d[[1]], d[[2]]), merge_dfs(d[[3]],d[[4]])), merge_dfs(d[[5]], d[[6]]))
  for (i in c(2,3,4,5,7)){precision_data[,i] <- as.numeric(precision_data[,i])}
    
  # Plot precision data
  # Compute positions for dashed lines
  x_start <- x_end <- y_start <- c()
  i_ <- 0.5
  for (i in 1:3){
    x_start <- c(x_start, rep(i_,2))
    x_end <- c(x_end, rep(i_+2,2))
    intra_study <- precision_data$intra.study[i*2]
    y_start <- c(y_start, rep(intra_study,2))
    i_ <- i_ + 2
  }
  y_end <- y_start
  
  # Plot
  fig6c_1 <- precision_data %>% 
    ggplot(mapping = aes(x = setting,
                         y = precision,
                         fill = dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1, vjust = 1),
          legend.direction = "horizontal") +
    ylim(c(0,1)) +
    scale_fill_manual(values = c("#332EBF", "#7D2EBF", "#BF2EAC")) +
    xlab("Setting") + ylab("Precision") + labs(fill = "Dataset") +
    geom_segment(aes(x = 0.5, xend = 6.5,
                     y = median(precision_data$precision), yend = median(precision_data$precision)),
                 linetype = "dashed", color = "#7a2eb9")
  
  fig6c_2 <- precision_data %>%
    ggplot(mapping = aes(x = setting,
                         y = placebo.swapped,
                         fill = dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1, vjust = 1),
          legend.direction = "vertical") +
    scale_fill_manual(values = c("#332EBF", "#7D2EBF", "#BF2EAC")) +
    xlab("Setting") + ylab("DE targets") + labs(fill = "Dataset") +
    geom_segment(aes(x = x_start, xend = x_end,
                     y = y_start,
                     yend = y_end),
                 linetype = "dashed", color = c("#332EBF","#332EBF", "#7D2EBF","#7D2EBF", "#BF2EAC","#BF2EAC"))

  legend <- get_legend(fig6c_1)
  g <- ggarrange(fig6c_1,fig6c_2,
                 ncol = 2,
                 common.legend = TRUE,
                 legend.grob = legend,
                 hjust = -2)
  title <- "Evaluation of swapping feasibility"
  fig6c <- annotate_figure(g, top = text_grob(title, vjust = 1, face = "bold"))
  
  
  # Plot
  fig6d_1 <- precision_data %>% 
    filter(setting %in% c("1 vs. 2", "1 vs. 3", "2 vs. 3")) %>% 
    ggplot(mapping = aes(x = setting,
                         y = top100,
                         fill = dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1, vjust = 1),
          legend.direction = "horizontal") +
    ylim(c(0,100)) +
    scale_fill_manual(values = c("#332EBF", "#7D2EBF", "#BF2EAC")) +
    xlab("Setting") + ylab("Overlap (top 100 genes)") + labs(fill = "Dataset") +
    geom_segment(aes(x = 0.5, xend = 3.5,
                     y = median(precision_data$top100), yend = median(precision_data$top100)),
                 linetype = "dashed", color = "#7a2eb9")
  
  fig6d_2 <- precision_data %>%
    ggplot(mapping = aes(x = setting,
                         y = overlap,
                         fill = dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1, vjust = 1),
          legend.direction = "horizontal") +
    ylim(c(0,1)) +
    scale_fill_manual(values = c("#332EBF", "#7D2EBF", "#BF2EAC")) +
    xlab("Setting") + ylab("Total overlap") + labs(fill = "Dataset") +
    geom_segment(aes(x = 0.5, xend = 6.5,
                     y = median(precision_data$overlap), yend = median(precision_data$overlap)),
                 linetype = "dashed", color = "#7a2eb9")
  
  legend <- get_legend(fig6d_2)
  g <- ggarrange(fig6d_1,fig6d_2,
                 ncol = 2,
                 common.legend = TRUE,
                 legend.grob = legend,
                 hjust = -2)
  title <- "Inter-study overlap"
  fig6d <- annotate_figure(g, top = text_grob(title, vjust = 1, face = "bold"))
  
  # will probably change this figure to be part mouse transcriptomics and part human proteomics
    # as I do not feel that the conclusions are strong enough for a full figure
  
  ## Put it all together
  fig6 <- ggarrange(fig6a,fig6b,fig6c,fig6d,
                    nrow = 4, heights = c(1.5,1.5,1,1),
                    labels = c("A", "B", "C", "D"))
  
  ## Save it
  fact <- 1.5
  ggsave("figures/fig6.png", plot = fig6, width = fact*(210-20), height = fact*(297-30),
         units = "mm", dpi = 185, bg = "white", limitsize = FALSE)
}





################################################################################
############################ Supplementary figures #############################
################################################################################

################################################################################
## Supplementary figure 15: Transcriptomics (intra-study volcano)
for (generate_figure in "sup_fig_15"){
  # Generate intra-study volcano plots
  myplots <- list(volcano_plot(short_title = TRUE),
                  volcano_plot(df1 = 2, short_title = TRUE),
                  volcano_plot(df1 = 3, short_title = TRUE))
  g <- ggarrange(plotlist = myplots,
                 nrow = length(myplots),
                 hjust = -2,
                 labels = c("A", "B", "C"))
  sup_fig_15 <- annotate_figure(g, top = text_grob("Intra-study volcano plots",
                                              vjust = 1.5, face = "bold"))
  
  fact <- 1.5
  ggsave("figures/supfig15.png", plot = sup_fig_15, width = fact*(210-20), height = fact*(297-30-35),
         units = "mm", dpi = 185, bg = "white")
}

################################################################################
## Supplementary figure 16: Transcriptomics (inter-study, chow)
for (generate_figure in "sup_fig_16"){
  # Generate inter-study volcano plots (lean chow)
  myplots <- list(volcano_plot(df2 = 2, type = "inter_chow", short_title = TRUE),
                  volcano_plot(df2 = 3, type = "inter_chow", short_title = TRUE),
                  volcano_plot(df1 = 2, df2 = 3, type = "inter_chow", short_title = TRUE))
  g <- ggarrange(plotlist = myplots,
                 nrow = length(myplots),
                 hjust = -2,
                 labels = c("A", "B", "C"))
  sup_fig_16 <- annotate_figure(g, top = text_grob("Inter-study volcano plots",
                                                   vjust = 1.5, face = "bold"))
  
  fact <- 1.5
  ggsave("figures/supfig16.png", plot = sup_fig_16, width = fact*(210-20), height = fact*(297-30-35),
         units = "mm", dpi = 185, bg = "white")
}

################################################################################
## Supplementary figure 17: Transcriptomics (inter-study, NASH)
for (generate_figure in "sup_fig_17"){
  # Generate inter-study volcano plots (lean chow)
  myplots <- list(volcano_plot(df2 = 2, type = "inter_nash", short_title = TRUE),
                  volcano_plot(df2 = 3, type = "inter_nash", short_title = TRUE),
                  volcano_plot(df1 = 2, df2 = 3, type = "inter_nash", short_title = TRUE))
  g <- ggarrange(plotlist = myplots,
                 nrow = length(myplots),
                 hjust = -2,
                 labels = c("A", "B", "C"))
  sup_fig_17 <- annotate_figure(g, top = text_grob("Inter-study volcano plots",
                                                   vjust = 1.5, face = "bold"))
  
  fact <- 1.5
  ggsave("figures/supfig17.png", plot = sup_fig_17, width = fact*(210-20), height = fact*(297-30-35),
         units = "mm", dpi = 185, bg = "white")
}

################################################################################
## Supplementary figure 18: Transcriptomics (placebo swaps)
for (generate_figure in "sup_fig_18"){
  # Generate placebo-swapped volcano plots
  myplots <- list(volcano_plot(df1 = 1, df2 = 2, type = "swap", short_title = TRUE),
                  volcano_plot(df1 = 1, df2 = 3, type = "swap", short_title = TRUE),
                  volcano_plot(df1 = 2, df2 = 1, type = "swap", short_title = TRUE),
                  volcano_plot(df1 = 2, df2 = 3, type = "swap", short_title = TRUE),
                  volcano_plot(df1 = 3, df2 = 1, type = "swap", short_title = TRUE),
                  volcano_plot(df1 = 3, df2 = 2, type = "swap", short_title = TRUE))
  g <- ggarrange(plotlist = myplots,
                 nrow = 3, ncol = 2,
                 hjust = -2,
                 labels = c("A", "B", "C", "D", "E", "F"))
  sup_fig_18 <- annotate_figure(g, top = text_grob("Placebo-swapped volcano plots",
                                                   vjust = 0.5, face = "bold"))
  
  fact <- 1.5
  ggsave("figures/supfig18.png", plot = sup_fig_18, width = fact*(210-20), height = fact*(297-30-40),
         units = "mm", dpi = 185, bg = "white")
}

################################################################################
## Supplementary figure 19: Transcriptomics (precision, overlap, etc.)




save("dge_t1_intra", "dge_t2_intra", "dge_t3_intra",
     "dge_t1_inter_chow_2", "dge_t1_inter_chow_3", "dge_t2_inter_chow_1",
     "dge_t2_inter_chow_3", "dge_t3_inter_chow_1", "dge_t3_inter_chow_2",
     "dge_t1_inter_nash_2", "dge_t1_inter_nash_3", "dge_t2_inter_nash_1",
     "dge_t2_inter_nash_3", "dge_t3_inter_nash_1", "dge_t3_inter_nash_2",
     "dge_t1_swap_2", "dge_t1_swap_3", "dge_t2_swap_1",
     "dge_t2_swap_3", "dge_t3_swap_1", "dge_t3_swap_2",
     "dge_intra", "dge_inter_chow", "dge_inter_nash", "dge_swap",
     file = "~/JDJG/data/curated/trans_dge.RData")
