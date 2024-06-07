### Initialize document
setwd("~/JDJG/clean/")
source("scripts/SomaScan_Mouse/adat_dens.R")
source("scripts/SomaScan_Mouse/baseline_analysis.R")
source("scripts/SomaScan_Mouse/pca.R")
source("scripts/SomaScan_Mouse/fc_plot.R")
source("scripts/SomaScan_Mouse/z_heatmap.R")
source("scripts/SomaScan_Mouse/volcano_function.R")
source("scripts/SomaScan_Mouse/precision_plot_function.R")
source("scripts/SomaScan_Mouse/bland_altman_function.R")
source("scripts/SomaScan_Mouse/bootstrap_plot_function.R")
source("scripts/SomaScan_Mouse/signif_targets_function.R")
source("scripts/SomaScan_Mouse/annotation_compass_function.R")
library(gridExtra); library(png); library(gridBase)

## Load data
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/bootstrap_fig7.RData")
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/bootstrap_supfig14.RData")
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/bootstrap_supfig16.RData")

################################################################################
################################################################################
################################# MAIN FIGURES #################################
################################################################################
################################################################################


################################################################################
############### FIGURE 3: Density plot, PCA, CV, and FC:FC plot ################
################################################################################
for (generate_figure in "FIGURE_3"){
  fact <- 2.5
  png(filename = "figures/fig3.png", width = fact*(210-20), height = fact*(297-30), units = "mm",
      pointsize = fact*7, res = 100)
  layout(mat = matrix(c(1,2,3,4,5,6,rep(7,4*3), rep(c(8,9,10),2)), nrow = 8, ncol = 3, byrow=TRUE))
  bot <- 2; top <- 0; title_pos <- 0.5
  # Density plots
  for (i in c(1,2)){
    for (norm in c("unnorm", "norm1", "norm2")){
      density_plot(i = i, norm = norm, bot = bot, top = top, title_pos = title_pos)
    }
  }
  # PCA and CV plots (combining base R and ggplot2)
  plot.new()
  vps <- baseViewports()
  pushViewport(vps$figure) ## think this moves us to the location of the plot
  vp1 <- plotViewport(c(1.8,1,0,1)) ## generate_figure new vp with margins (play with the values)
  print(pca_cv_final, vp = vp1)
  # FC plots
  fc_plot(handle = "init_",
          text = "No outliers removed", height = 800, width = 1067,
          top100 = FALSE, save = FALSE)
  fc_plot(handle = "init_1_5_",
          text = "No outliers removed, medNorm", height = 800, width = 1067,
          top100 = FALSE, save = FALSE)
  fc_plot(handle = "log_",
          text = "Outliers removed", height = 800, width = 1067,
          top100 = FALSE, save = FALSE)
  dev.off()
}



################################################################################
############## FIGURE 4: Baseline volcano plots + precision plot ###############
################################################################################
for (generate_figure in "FIGURE_4"){
  
  ## Volcano plots
  # Create intra-study case:control plots
  myplots <- list(volcano_plot(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               title = "SomaMouse1", case = 3),
                  volcano_plot(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                               title = "SomaMouse2", case = 2))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc1 <- annotate_figure(g, top = text_grob("DIO-NASH Vehicle vs. Chow Vehicle (intra-study)",
                                              vjust = 1.5, face = "bold"))
  
  # Create inter-study control:control plots
  myplots <- list(volcano_plot(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 3, group2 = 2, title = "DIO-NASH Vehicle"),
                  volcano_plot(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 1, group2 = 1, title = "Lean chow Vehicle"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc2 <- annotate_figure(g, top = text_grob("SomaMouse1 vs. SomaMouse2 (control:control)",
                                              vjust = 1.5, face = "bold"))
  
  # Create placebo-swapped plots
  myplots <- list(volcano_plot(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 3, group2 = 1, title = "SomaMouse1 (case) vs. SomaMouse2 (control)"),
                  volcano_plot(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                               soma2 = log_JYNR210602, meta2 = meta_mouse_JYNR210602,
                               group1 = 2, group2 = 1, title = "SomaMouse2 (case) vs. SomaMouse1 (control)"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc3 <- annotate_figure(g, top = text_grob("DIO-NASH Vehicle vs. Chow Vehicle (swapped)",
                                              vjust = 1.5, face = "bold"))
  
  
  ## Precision plot
  pre_plot <- precision_plot(soma1_list = list(log_JYNR210602),
                             meta1_list = list(meta_mouse_JYNR210602),
                             soma2_list = list(log_JYNR220601),
                             meta2_list = list(meta_mouse_JYNR220601),
                             norm_methods = c("Placebo swapped"),
                             x_angle = 0, x_hjust = 0.5)
  
  ## Put it all together
  fig4 <- ggarrange(volc1,NULL,volc2,NULL,volc3,NULL,pre_plot,
                    nrow = 7, heights = c(1,0.05,1,0.05,1,0.05,1),
                    labels = c("A", "", "B", "", "C", "", "D"))
  
  
  ## Save it
  fact <- 1.5
  ggsave("figures/fig4.png", plot = fig4, width = fact*(210-20), height = fact*(297-30),
         units = "mm", dpi = 185, bg = "white")
}





################################################################################
#### FIGURE 5: Cross-norm and limma: Bland-Altman, volcano, PCA, precision #####
################################################################################
for (generate_figure in "FIGURE_5"){
  
  ## Bland-Altman plots
  ba_plots <- list()
  norms <- c(FALSE, TRUE, "limma")
  comps <- c("case", "control")
  for (i in 1:2){
    comp <- comps[i]
    myplots <- list()
    limits <- determine_limits(soma1_dfs = list(log_JYNR210602, norm_JYNR210602, limma_JYNR210602),
                               soma2_dfs = list(log_JYNR220601, norm_JYNR220601, limma_JYNR220601),
                               comp = comp)
    xlim <- c(limits[1], limits[2]); ylim <- c(limits[3], limits[4])
    for (norm in norms){
      # Define dataframe
      if (norm == TRUE){
        soma1 <- "crossnorm_JYNR210602"; soma2 <- "crossnorm_JYNR220601"
        norm <- "Cross-normalised (SomaLogic)"
        print("cross-norm")
      } else if (norm == "limma"){
        soma1 <- "limma_JYNR210602"; soma2 <- "limma_JYNR220601"
        norm <- "Batch effects corrected (limma)"
        print ("limma")
      } else {
        soma1 <- "log_JYNR210602"; soma2 <- "log_JYNR220601"
        norm <- "Baseline"
        print("log")
      }
      soma1 <- get(paste0(soma1)); soma2 <- get(paste0(soma2))
      
      # Plot
      plt <- list(bland_altman(soma1 = soma1, soma2 = soma2, normalisation = norm,
                               comparison = comp, title = norm,
                               xlim = xlim, ylim = ylim, short_y_lab = TRUE))
      myplots <- append(myplots, plt)
    }
    g <- ggarrange(plotlist = myplots,
                   ncol = length(norms),
                   nrow = 1,
                   hjust = -2)
    title <- paste0(str_to_title(comp), " vs. ", comp)
    ba_plot <- annotate_figure(g, top = text_grob(title, vjust = 1, face = "italic"))
    ba_plots <- append(ba_plots, list(ba_plot))
  }
  ba_final <- annotate_figure(ggarrange(ba_plots[[1]], NULL, ba_plots[[2]],
                                        nrow = 3, ncol = 1, heights = c(1,0.05,1)),
                              top = text_grob("Bland-Altman plot, SomaMouse1 vs. SomaMouse2",
                                              vjust = 1, face = "bold"))
  
  
  ## Volcano plots
  # Create inter-study control:control plots for cross-norm and limma
  legend <- get_legend(volcano_plot(soma1 = crossnorm_JYNR210602, meta1 = meta_mouse_JYNR210602,
                                    soma2 = crossnorm_JYNR220601, meta2 = meta_mouse_JYNR220601,
                                    group1 = 1, group2 = 1, title = "Cross-normalisation (SomaLogic)"))
  myplots <- list(volcano_plot(soma1 = crossnorm_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = crossnorm_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 1, group2 = 1, title = "Cross-normalisation (SomaLogic)"),
                  volcano_plot(soma1 = limma_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = limma_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 1, group2 = 1, title = "Batch effect correction (limma)",
                               title_vjust = 3))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2,
                 common.legend = TRUE,
                 legend.grob = legend,
                 legend = "top")
  volc <- annotate_figure(g, top = text_grob("SomaMouse1 vs. SomaMouse2: Lean chow vehicle",
                                             vjust = 1, face = "bold"))
  
  
  ## Principal component analysis
  myplots <- list(pcaplot(crossnorm_data, title = "Cross-normalisation (SomaLogic)",
                          legend = "none", colour = "Treatment", point_size = 3,
                          short_labs = TRUE),
                  pcaplot(limma_data, title = "Batch effect correction (limma)",
                          legend = "none", colour = "Treatment", point_size = 3,
                          short_labs = TRUE))
  legend <- pcaplot(baseline_data, colour = "Treatment", legend = "top", get_legend = TRUE)
  pca <- ggarrange(plotlist = myplots,
                   nrow = length(myplots),
                   common.legend = TRUE,
                   legend.grob = legend,
                   legend = "left",
                   hjust = -2)
  
  ## Precision plot
  pre_plot <- precision_plot(soma1_list = list(log_JYNR210602, crossnorm_JYNR210602, limma_JYNR210602),
                             meta1_list = list(meta_mouse_JYNR210602, meta_mouse_JYNR210602, meta_mouse_JYNR210602),
                             soma2_list = list(log_JYNR220601, crossnorm_JYNR220601, limma_JYNR220601),
                             meta2_list = list(meta_mouse_JYNR220601, meta_mouse_JYNR220601, meta_mouse_JYNR220601),
                             norm_methods = c("Baseline", "Crossnorm", "limma"),
                             short_y_lab = TRUE, fixed_ref = TRUE,
                             x_angle = 0, x_hjust = 0.5)
  
  ## Put it all together
  fig5 <- ggarrange(ba_final,NULL,volc,NULL,pca,NULL,pre_plot,
                    nrow = 7, heights = c(1,0.05,0.5,0.05,1,0.05,0.5),
                    labels = c("A", "", "B", "", "C", "", "D"))
  
  
  ## Save it
  fact <- 1.5
  ggsave("figures/fig5.png", plot = fig5, width = fact*(210-20), height = fact*(297-30),
         units = "mm", dpi = 185, bg = "white")
}



################################################################################
######## FIGURE 6: Realistic scenarios, calibrator-based normalisation #########
################################################################################
for (generate_figure in "FIGURE_6"){
  
  ## Schematic of realistic scenarios
  legend1 <- pcaplot(list(log_JYNR210602, comb_meta[comb_meta$Dataset == "SomaMouse1",]),
                     colour = "Treatment", legend = "top", get_legend = TRUE)
  legend2 <- pcaplot(list(log_JYNR220601, comb_meta[comb_meta$Dataset == "SomaMouse2",]),
                     colour = "Treatment", legend = "top", get_legend = TRUE)
  legend3 <- pcaplot(baseline_data, colour = "Treatment", get_legend = TRUE)
  ggarrange(legend1, legend2, legend3, ncol = 3)
  ggsave("figures/fig6a_scrambled.png", width = 200, height = 200, units = "mm")
  
  ## Calibrator-based normalisation
  # Create inter-study control:control plots (baseline)
  myplots <- list(volcano_plot(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 3, group2 = 2, title = "DIO-NASH Vehicle"),
                  volcano_plot(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 1, group2 = 1, title = "Lean chow Vehicle"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc_cc_base <- annotate_figure(g, top = text_grob("Baseline",
                                                     vjust = 1.5, face = "italic"))
  
  # Generate inter-study control:control plots (calibrator-based normalisation)
  myplots <- list(volcano_plot(soma1 = cal_indi_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 3, group2 = 2, title = "DIO-NASH Vehicle"),
                  volcano_plot(soma1 = cal_indi_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 1, group2 = 1, title = "Lean chow Vehicle"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc_cc_cal <- annotate_figure(g, top = text_grob("Calibrator-based cross-normalisation",
                                                    vjust = 1.5, face = "italic"))
  g <- ggarrange(volc_cc_base, volc_cc_cal, ncol = 1, nrow = 2)
  volc_cc_final <- annotate_figure(g, top = text_grob("Control:control volcano plots",
                                                      vjust = 1, face = "bold"))
  
  ## Precision plot
  pre_plot <- precision_plot(
    cases = list(
      list(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # No normalisation (1-4)
           crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Cross-normalisation
           case_control_norm_dfs[[1]][[1]], # scenario 1, cross-normalised (JYNR210602)
           case_control_norm_dfs[[5]][[1]], # scenario 2, cross-normalised (JYNR210602)
           limma_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Limma (1-4)
           case_control_norm_dfs[[3]][[1]], # scenario 1, limma (JYNR210602)
           case_control_norm_dfs[[7]][[1]], # scenario 2, limma (JYNR210602)
           cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 3,]), # Calibrator-normalised (indi SFs)
      list(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
           crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Cross-normalisation
           case_control_norm_dfs[[2]][[1]], # scenario 1, cross-normalised (JYNR220601)
           case_control_norm_dfs[[6]][[1]], # scenario 2, cross-normalised (JYNR220601)
           limma_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Limma (1-4)
           case_control_norm_dfs[[4]][[1]], # scenario 1, limma (JYNR220601)
           case_control_norm_dfs[[8]][[1]], # scenario 2, limma (JYNR220601)
           log_JYNR220601[meta_mouse_JYNR220601$Group == 2,])), #  Calibrator-normalised (indi SFs)
    placebos = list(
      list(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
           crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Cross-normalisation
           case_control_norm_dfs[[1]][[2]], # scenario 1, cross-normalised (JYNR210602)
           case_control_norm_dfs[[5]][[2]], # scenario 2, cross-normalised (JYNR210602)
           limma_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Limma (1-4)
           case_control_norm_dfs[[3]][[2]], # scenario 1, limma (JYNR210602)
           case_control_norm_dfs[[7]][[2]], # scenario 2, limma (JYNR210602)
           log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]), #  Calibrator-normalised (indi SFs)
      list(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # No normalisation (1-4)
           crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Cross-normalisation
           case_control_norm_dfs[[2]][[2]], # scenario 1, cross-normalised (JYNR220601)
           case_control_norm_dfs[[6]][[2]], # scenario 2, cross-normalised (JYNR220601)
           limma_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Limma (1-4)
           case_control_norm_dfs[[4]][[2]], # scenario 2, limma (JYNR220601)
           case_control_norm_dfs[[8]][[2]], # scenario 2, limma (JYNR220601)
           cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 1,])), # Calibrator-normalised (indi SFs)
    norm_methods = c("Baseline", "Cross", "Cross, real 1", "Cross, real 2",
                     "limma", "limma, real 1", "limma, real 2", "Cal"),
    short_y_lab = FALSE, fixed_ref = TRUE, #just_numbers = TRUE,
    x_angle = 30, x_hjust = 1, x_vjust = 1, meta_aligned = TRUE)
  
  
  ## Put it all together
  fig6_bc <- ggarrange(volc_cc_final,NULL,pre_plot,
                       nrow = 3, heights = c(1,0.05,0.5),
                       labels = c("B", "", "C"))
  
  
  ## Save it
  fact <- 1.5
  ggsave("figures/fig6bc.png", plot = fig6_bc, width = fact*(210-20), height = fact*(181),
         units = "mm", dpi = 185, bg = "white")
}


################################################################################
########################## FIGURE 7: Bootstrapping #############################
################################################################################
for (generate_figure in "FIGURE_7"){
  
  ## Fig 7a
  fig7a <- bootstrap_plot(bootstrap_stats = bootstrap_stats_fig7a, n1 = c(15,10), n2 = c(7,5))
  
  ## Fig 7b
  fig7b_data <- bootstrap_n_test_fig7b_100
  
  # Plot it
  fig7b <- fig7b_data %>% 
    ggplot(mapping = aes(x = dataset,
                         y = sig_tar,
                         fill = bootstrap,
                         label = n_c)) +
    geom_boxplot()+
    ylab("DE targets (intra-study)") +
    # ylab("DE targets") +
    xlab("Dataset") +
    ggtitle(paste0("Bootstrapping test for different sampling sizes")) +
    theme(plot.title = element_text(face = "bold")) +
    annotate("text", label = paste0("it = ", unique(fig7b_data$it)),
             x = 0.5, y = max(c(fig7b_data$ci_up, fig7b_data$sig_tar)), hjust=0, vjust=1) +
    labs(fill = "Bootstrap") +
    scale_fill_manual(values = c("#FFA355", "#C6A88E", "#C6A88E", "#C6A88E", "#C6A88E",
                                 "#99ABBB", "#99ABBB", "#99ABBB", "#55B1FF")) +
    geom_hline(yintercept = 479, linetype = "dashed", col = "#FFA355") +
    geom_hline(yintercept = 1924, linetype = "dashed", col = "#55B1FF") +
    annotate(geom="text", x = 1, y = median(fig7b_data$sig_tar[201:300]), label = paste0())
  
  # Fig 7c
  fig7c <- bootstrap_plot(bootstrap_stats = bootstrap_stats_fig7c, n1 = c(45,30), n2 = c(4,3))
  
  ## Put it all together
  fig7 <- ggarrange(fig7a,NULL,fig7b,NULL,fig7c,
                    nrow = 5, heights = c(1,0.05,1,0.05,1),
                    labels = c("A", "", "B", "", "C"))
  
  ## Save it
  fact <- 1.5
  ggsave("figures/fig7.png", plot = fig7, width = fact*(210-20), height = fact*(297-30),
         units = "mm", dpi = 185, bg = "white")
}
# median(c(bootstrap_stats_fig7a[[1]]$precision.1, bootstrap_stats_fig7a[[1]]$precision.2)) # 56.75%
# median(c(bootstrap_stats_fig7d[[1]]$precision.1, bootstrap_stats_fig7d[[1]]$precision.2)) # 56.89%







################################################################################
################################################################################
### Supplementary figures
## Supplementary figure 1: Density plots (all outliers vs. 5 outliers removed)
for (generate_figure in "sup_fig_1"){
  fact <- 2.5 * 2
  png(filename = "figures/supfig1.png", width = fact*(210-20), height = fact*(297-30-50), units = "mm",
      pointsize = fact*7, res = 100)
  layout(mat = matrix(c(1:12), nrow = 4, ncol = 3, byrow=TRUE))
  bot <- 0; top <- 2; title_pos <- 0.5
  for (i in c(1,2)){
    if (i == 2){
      bot <- 2; top <- 0
    }
    for (norm in c("unnorm", "norm1", "norm2")){
      
      density_plot(i = i, norm = norm, bot = bot, top = top, title_pos = title_pos)
    }
  }
  bot <- 0; top <- 2
  for (i in c(1,2)){
    if (i == 2){
      bot <- 2; top <- 0
    }
    for (norm in c("unnorm_outrm", "norm1_outrm", "norm2_outrm")){
      density_plot(i = i, norm = norm, bot = bot, top = top, title_pos = title_pos)
    }
  }
  dev.off()
}


################################################################################
## Supplementary figure 2: PCA
for (generate_figure in "sup_fig_2"){
  pcaplot(log_1_5_data, title = "PCA: Outliers removed, interplate median signal normalisation", legend = "left", colour = "Treatment")
  ggsave("figures/supfig2.png", width = 2800, height = 500, units = "px", dpi = 150, bg = "white")
  
}


################################################################################
## Supplementary figure 3: CV plots
for (generate_figure in "sup_fig_3"){
  myplots <- list(cv_plot(init_JYNR210602, init_JYNR220601, init_meta_JYNR210602, init_meta_JYNR220601,
                          text = "No outliers removed",x_max = 15, hide_legend = TRUE),
                  cv_plot(log_JYNR210602, log_JYNR220601, meta_mouse_JYNR210602, meta_mouse_JYNR220601,
                          text = "Outliers removed", x_max = 15, hide_legend = TRUE),
                  cv_plot(init_1_5_JYNR210602, init_1_5_JYNR220601, init_meta_JYNR210602, init_meta_JYNR220601,
                          text = "No outliers removed, medNorm", x_max = 15, hide_legend = TRUE),
                  cv_plot(log_1_5_JYNR210602, log_1_5_JYNR220601, meta_mouse_JYNR210602, meta_mouse_JYNR220601,
                          text = "Outliers removed, medNorm", x_max = 15, hide_legend = TRUE)
  )
  legend <- get_legend(cv_plot(init_JYNR210602, init_JYNR220601, init_meta_JYNR210602,
                               init_meta_JYNR220601, only_cont = TRUE))
  cv_final <- ggarrange(plotlist = myplots,
                        nrow = 2, ncol = 2,
                        common.legend = TRUE,
                        legend.grob = legend,
                        legend = "bottom",
                        hjust = -2)
  cv_final
  fact <- 2.5
  ggsave("figures/supfig3.png", width = fact*(210-20), height = fact*(297-30-150),
         units = "mm", dpi = 500, bg = "white")
}


################################################################################
## Supplementary figure 4: FC:FC plot
for (generate_figure in "sup_fig_4"){
  fact <- 1
  scale <- 0.5
  png(filename = "figures/supfig4.png", width = fact*(210-20), height = fact*(297-30-50), units = "mm",
      pointsize = fact*7*1.5, res = 150)
  layout(mat = matrix(c(1:4), nrow = 2, ncol = 2, byrow=TRUE))
  fc_plot(handle = "init_",
          text = "No outliers removed", height = 800*scale, width = 1067*scale,
          top100 = FALSE, save = FALSE)
  fc_plot(handle = "log_",
          text = "Outliers removed", height = 800, width = 1067,
          top100 = FALSE, save = FALSE)
  fc_plot(handle = "init_1_5_",
          text = "No outliers removed, medNorm", height = 800, width = 1067,
          top100 = FALSE, save = FALSE)
  fc_plot(handle = "log_1_5_",
          text = "Outliers removed, medNorm", height = 800, width = 1067,
          top100 = FALSE, save = FALSE)
  dev.off()
}


################################################################################
## Supplementary figure 5: Heatmaps
for (generate_figure in "sup_fig_5"){
  handles <- c("init_", "log_", "init_1_5_", "log_1_5_", "z_joint_")
  subfig <- c("a", "b", "c", "d", "e")
  legends <- c(rep(FALSE,4), TRUE)
  width <- c(rep(12,4),14); height <- rep(9,5)
  for (i in 1:5){
    handle <- handles[i]
    soma_data_1 <- t(as.matrix(get(paste0(handle, "JYNR210602"))))
    soma_data_2 <- t(as.matrix(get(paste0(handle, "JYNR220601"))))
    soma_data_merged <- merge(soma_data_1, soma_data_2, by = "row.names", all = TRUE)
    rownames(soma_data_merged) <- soma_data_merged[,1]
    soma_data_merged <- soma_data_merged[-1]
    
    # Get metadata
    if (handle %in% c("init_", "init_1_5_")){
      meta1 <- init_meta_JYNR210602; meta2 <- init_meta_JYNR220601
    } else if (handle %in% c("log_", "log_1_5_", "z_joint_")){
      meta1 <- meta_mouse_JYNR210602; meta2 <- meta_mouse_JYNR220601
    }
    
    meta1_cc <- meta1[meta1$Group %in% c(1,3), colnames(meta1) %in% colnames(meta2)]
    meta2_cc <- meta2[meta2$Group %in% c(1,2), colnames(meta2) %in% colnames(meta1)]
    grps <- c(rownames(meta1[meta1$Group %in% c(1,3),]),
              rownames(meta2[meta2$Group %in% c(1,2),]))
    # Keep only case/control vehicle groups
    soma_data_merged <- soma_data_merged[,colnames(soma_data_merged) %in% grps]
    meta_merged <- dplyr::full_join(meta1_cc, meta2_cc)
    colnames(meta_merged)[colnames(meta_merged) == "dataset"] <- "Dataset"
    # meta_merged$Group[meta_merged$id == "118499"] <- 9 # comment out
    
    sampleinfo <- meta_merged %>%
      select(Model, Dataset)
    corMatrix <- cor(soma_data_merged,use = "c")
    rownames(sampleinfo) <- colnames(corMatrix)
    
    
    ann_colours <- list(Dataset = c(SomaMouse1 = "#FFA355",
                                    SomaMouse2 = "#55B1FF"),
                        Model = c(Healthy = "#53BF2E", NASH = "#BF3D35"))
    if (handle == "log_"){
      title <- paste0("Outliers removed")
    }
    if (handle == "init_"){
      title <- paste0("No outliers removed")
    }
    if (handle == "log_1_5_"){
      title <- paste0("Outliers removed, median normalised")
    }
    if (handle == "init_1_5_"){
      title <- paste0("No outliers removed, median normalised")
    }
    if (handle == "z_joint_"){
      title <- paste0("Outliers removed, Z-scores")
    }
    filename <- paste0("figures/supfig5", subfig[i], ".png")
    fact <- 0.9
    pheatmap(corMatrix, annotation_col = sampleinfo,
             show_rownames = FALSE, show_colnames = FALSE, annotation_colors = ann_colours,
             filename = filename, legend = FALSE, annotation_legend = legends[i],
             main = title, width = width[i]*fact, height = height[i]*fact, fontsize = 16)
  }
}
# put them together to one final figure manually


################################################################################
## Supplementary figure 6: Overlap (top100)
for (generate_figure in "sup_fig_6"){
  overlaps <- list(top100_overlap(init_JYNR210602, init_JYNR220601, init_meta_JYNR210602, init_meta_JYNR220601),
                   top100_overlap(log_JYNR210602, log_JYNR220601, meta_mouse_JYNR210602, meta_mouse_JYNR220601),
                   top100_overlap(init_1_5_JYNR210602, init_1_5_JYNR220601, init_meta_JYNR210602, init_meta_JYNR220601),
                   top100_overlap(log_1_5_JYNR210602, log_1_5_JYNR220601, meta_mouse_JYNR210602, meta_mouse_JYNR220601))
  norm_methods <- c("Outliers present", "Outliers removed", "Outliers present, medNorm", "Outliers removed, medNorm")
  overlap_df <- tibble(`Normalisation method` = factor(norm_methods, levels = norm_methods),
                       Overlap = unlist(overlaps))
  p1 <- overlap_df %>% 
    ggplot(mapping = aes(x = `Normalisation method`,
                         y = Overlap,
                         fill = `Normalisation method`)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 30,
                                     hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("Overlap between datasets,\ntop 100 significant targets by absolute fold change") +
    ylab("Overlap (top 100 targets)") +
    ylim(c(0,100))
  sup_fig_6 <- ggarrange(p1, legend = "none")
  ggsave("figures/supfig6.png", plot = sup_fig_6, width = fact*(210-20), height = fact*(297-30-160),
         units = "mm", dpi = 500, bg = "white")
}


################################################################################
## Supplementary figure 7: Precision
## Precision plot
for (generate_figure in "sup_fig_7"){
  sup_fig_7 <- precision_plot(soma1_list = list(init_JYNR210602, log_JYNR210602, init_1_5_JYNR210602, log_1_5_JYNR210602),
                              meta1_list = list(init_meta_JYNR210602, meta_mouse_JYNR210602, init_meta_JYNR210602, meta_mouse_JYNR210602),
                              soma2_list = list(init_JYNR220601, log_JYNR220601, init_1_5_JYNR220601, log_1_5_JYNR220601),
                              meta2_list = list(init_meta_JYNR220601, meta_mouse_JYNR220601, init_meta_JYNR220601, meta_mouse_JYNR220601),
                              norm_methods = c("No outliers removed",
                                               "Outliers removed",
                                               "No outliers removed, medNorm",
                                               "Outliers removed, medNorm"))
  
  ## Save it
  fact <- 1.5
  ggsave("figures/supfig7.png", plot = sup_fig_7, width = fact*(210-20), height = fact*(297-30-160),
         units = "mm", dpi = 185, bg = "white")
}


################################################################################
## Supplementary figure 8: Volcano (intra-study)
for (generate_figure in "sup_fig_8"){
  ## Volcano plots
  # Generate_figure intra-study case:control plots (baseline)
  myplots <- list(volcano_plot(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               title = "SomaMouse1", case = 3),
                  volcano_plot(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                               title = "SomaMouse2", case = 2))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc1 <- annotate_figure(g, top = text_grob("Baseline",
                                              vjust = 1.5, face = "bold"))
  
  # Generate_figure intra-study case:control plots (baseline)
  myplots <- list(volcano_plot(soma1 = crossnorm_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               title = "SomaMouse1", case = 3),
                  volcano_plot(soma1 = crossnorm_JYNR220601, meta1 = meta_mouse_JYNR220601,
                               title = "SomaMouse2", case = 2))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc2 <- annotate_figure(g, top = text_grob("Cross-normalisation (SomaLogic)",
                                              vjust = 1.5, face = "bold"))
  
  # Generate_figure intra-study case:control plots (baseline)
  myplots <- list(volcano_plot(soma1 = limma_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               title = "SomaMouse1", case = 3),
                  volcano_plot(soma1 = limma_JYNR220601, meta1 = meta_mouse_JYNR220601,
                               title = "SomaMouse2", case = 2))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc3 <- annotate_figure(g, top = text_grob("Batch effect correction (limma)",
                                              vjust = 1.5, face = "bold"))
  
  ## Put it all together
  g <- ggarrange(volc1,NULL,volc2,NULL,volc3,
                 nrow = 5, heights = c(1,0.05,1,0.05,1),
                 labels = c("A", "", "B", "", "C", "", "D"))
  sup_fig_8 <- annotate_figure(g, top = text_grob("DIO-NASH Vehicle vs. Chow Vehicle (intra-study)",
                                                  face = "bold"))
  
  ## Save it
  fact <- 1.5
  ggsave("figures/supfig8.png", plot = sup_fig_8, width = fact*(210-20), height = fact*(297-30-50),
         units = "mm", dpi = 185, bg = "white")
}


################################################################################
## Supplementary figure 9: Volcano (control:control)
for (generate_figure in "sup_fig_9"){
  ## Volcano plots
  # Generate_figure inter-study control:control plots (baseline)
  myplots <- list(volcano_plot(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 3, group2 = 2, title = "DIO-NASH Vehicle"),
                  volcano_plot(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 1, group2 = 1, title = "Lean chow Vehicle"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc1 <- annotate_figure(g, top = text_grob("Baseline",
                                              vjust = 1.5, face = "bold"))
  
  # Generate_figure inter-study control:control plots (crossnorm)
  myplots <- list(volcano_plot(soma1 = crossnorm_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = crossnorm_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 3, group2 = 2, title = "DIO-NASH Vehicle"),
                  volcano_plot(soma1 = crossnorm_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = crossnorm_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 1, group2 = 1, title = "Lean chow Vehicle"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc2 <- annotate_figure(g, top = text_grob("Cross-normalisation (SomaLogic)",
                                              vjust = 1.5, face = "bold"))
  
  # Generate_figure inter-study control:control plots (limma)
  myplots <- list(volcano_plot(soma1 = limma_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = limma_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 3, group2 = 2, title = "DIO-NASH Vehicle"),
                  volcano_plot(soma1 = limma_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = limma_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 1, group2 = 1, title = "Lean chow Vehicle"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc3 <- annotate_figure(g, top = text_grob("Batch effect correction (limma)",
                                              vjust = 1.5, face = "bold"))
  
  ## Put it all together
  g <- ggarrange(volc1,NULL,volc2,NULL,volc3,
                 nrow = 5, heights = c(1,0.05,1,0.05,1),
                 labels = c("A", "", "B", "", "C", "", "D"))
  sup_fig_9 <- annotate_figure(g, top = text_grob("Control:control (inter-study)",
                                                  face = "bold"))
  
  ## Save it
  fact <- 1.5
  ggsave("figures/supfig9.png", plot = sup_fig_9, width = fact*(210-20), height = fact*(297-30-50),
         units = "mm", dpi = 185, bg = "white")
}


################################################################################
## Supplementary figure 10: Volcano (placebo swaps)
for (generate_figure in "sup_fig_10"){
  ## Volcano plots
  # Generate_figure placebo-swapped plots (baseline)
  myplots <- list(volcano_plot(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 3, group2 = 1, title = "SomaMouse1 (case) vs. SomaMouse2 (control)"),
                  volcano_plot(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                               soma2 = log_JYNR210602, meta2 = meta_mouse_JYNR210602,
                               group1 = 2, group2 = 1, title = "SomaMouse2 (case) vs. SomaMouse1 (control)"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc1 <- annotate_figure(g, top = text_grob("Baseline",
                                              vjust = 1.5, face = "bold"))
  
  # Generate_figure placebo-swapped plots (crossnorm)
  myplots <- list(volcano_plot(soma1 = crossnorm_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = crossnorm_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 3, group2 = 1, title = "SomaMouse1 (case) vs. SomaMouse2 (control)"),
                  volcano_plot(soma1 = crossnorm_JYNR220601, meta1 = meta_mouse_JYNR220601,
                               soma2 = crossnorm_JYNR210602, meta2 = meta_mouse_JYNR210602,
                               group1 = 2, group2 = 1, title = "SomaMouse2 (case) vs. SomaMouse1 (control)"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc2 <- annotate_figure(g, top = text_grob("Cross-normalisation (SomaLogic)",
                                              vjust = 1.5, face = "bold"))
  
  # Generate_figure placebo-swapped plots (limma)
  myplots <- list(volcano_plot(soma1 = limma_JYNR210602, meta1 = meta_mouse_JYNR210602,
                               soma2 = limma_JYNR220601, meta2 = meta_mouse_JYNR220601,
                               group1 = 3, group2 = 1, title = "SomaMouse1 (case) vs. SomaMouse2 (control)"),
                  volcano_plot(soma1 = limma_JYNR220601, meta1 = meta_mouse_JYNR220601,
                               soma2 = limma_JYNR210602, meta2 = meta_mouse_JYNR210602,
                               group1 = 2, group2 = 1, title = "SomaMouse2 (case) vs. SomaMouse1 (control)"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 hjust = -2)
  volc3 <- annotate_figure(g, top = text_grob("Batch effect correction (limma)",
                                              vjust = 1.5, face = "bold"))
  
  ## Put it all together
  g <- ggarrange(volc1,NULL,volc2,NULL,volc3,
                 nrow = 5, heights = c(1,0.05,1,0.05,1),
                 labels = c("A", "", "B", "", "C", "", "D"))
  sup_fig_10 <- annotate_figure(g, top = text_grob("DIO-NASH Vehicle vs. Chow Vehicle (swapped)",
                                                   face = "bold"))
  
  ## Save it
  fact <- 1.5
  ggsave("figures/supfig10.png", plot = sup_fig_10, width = fact*(210-20), height = fact*(297-30-50),
         units = "mm", dpi = 185, bg = "white")
}


################################################################################
## Supplementary figure 11: PCA (crossnorm + limma, only vehicle groups)
for (generate_figure in "sup_fig_11"){
  ## Principal component analysis
  myplots <- list(pcaplot(crossnorm_data, title = "Cross-normalisation (SomaLogic)",
                          legend = "none", colour = "Treatment", point_size = 3,
                          groups = c("Chow Vehicle", "DIO-NASH Vehicle"), short_labs = TRUE),
                  pcaplot(limma_data, title = "Batch effect correction (limma)",
                          legend = "none", colour = "Treatment", point_size = 3,
                          groups = c("Chow Vehicle", "DIO-NASH Vehicle"), short_labs = TRUE))
  legend <- pcaplot(baseline_data, colour = "Treatment", legend = "top", get_legend = TRUE,
                    groups = c("Chow Vehicle", "DIO-NASH Vehicle"))
  sup_fig_11 <- ggarrange(plotlist = myplots,
                          nrow = length(myplots),
                          common.legend = TRUE,
                          legend.grob = legend,
                          legend = "left",
                          hjust = -2)
  ## Save it
  fact <- 1.5
  ggsave("figures/supfig11.png", plot = sup_fig_11, width = fact*(210-20), height = fact*(297-30-180),
         units = "mm", dpi = 185, bg = "white")
}

################################################################################
## Supplementary figure 12: Volcano plots (realistic, limma, swaps)
for (generate_figure in "sup_fig_12"){
  # Generate_figure placebo-swapped plots (realistic scenario 1)
  myplots <- list(volcano_plot(soma1 = case_control_norm_dfs[[3]][[1]], # scenario 1, limma (JYNR210602, case)
                               soma2 = case_control_norm_dfs[[3]][[2]], # scenario 1, limma (JYNR220601, placebo)
                               meta_aligned = TRUE, title = "SomaMouse1 (case) vs. SomaMouse2 (control)"),
                  volcano_plot(soma1 = case_control_norm_dfs[[4]][[1]], # scenario 1, limma (JYNR220601, case)
                               soma2 = case_control_norm_dfs[[4]][[2]], # scenario 2, limma (JYNR210602, case)
                               meta_aligned = TRUE, title = "SomaMouse2 (case) vs. SomaMouse1 (control)"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 common.legend = TRUE,
                 hjust = -2)
  volc1 <- annotate_figure(g, top = text_grob("Scenario 1",
                                              vjust = 1.5, face = "bold"))
  
  # Generate_figure placebo-swapped plots (limma)
  myplots <- list(volcano_plot(soma1 = case_control_norm_dfs[[7]][[1]], # scenario 2, limma (JYNR210602, case)
                               soma2 = case_control_norm_dfs[[7]][[2]], # scenario 2, limma (JYNR220601, placebo)
                               meta_aligned = TRUE, title = "SomaMouse1 (case) vs. SomaMouse2 (control)"),
                  volcano_plot(soma1 = case_control_norm_dfs[[8]][[1]], # scenario 2, limma (JYNR220601, case)
                               soma2 = case_control_norm_dfs[[8]][[2]], # scenario 2, limma (JYNR210602, case)
                               meta_aligned = TRUE, title = "SomaMouse2 (case) vs. SomaMouse1 (control)"))
  g <- ggarrange(plotlist = myplots,
                 ncol = length(myplots),
                 common.legend = TRUE,
                 hjust = -2)
  volc2 <- annotate_figure(g, top = text_grob("Scenario 2",
                                              vjust = 1.5, face = "bold"))
  
  ## Put it all together
  g <- ggarrange(volc1,volc2,
                 nrow = 2, heights = c(0.6,1),
                 labels = c("A", "B"))
  sup_fig_12 <- annotate_figure(g, top = text_grob("Volcano plots: Realistic scenarios with limma",
                                                   face = "bold"))
  
  ## Save it
  fact <- 1.5
  ggsave("figures/supfig12.png", plot = sup_fig_12, width = fact*(210-20), height = fact*(297-30-140),
         units = "mm", dpi = 185, bg = "white")
}


################################################################################
## Supplementary figure 13: Aptamer quality filtering (precision plot)
for (generate_figure in "sup_fig_13"){
  ## Precision plots
  # Baseline
  pre_plot_b <- precision_plot(
    cases = list(
      list(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # No normalisation (1-4)
           crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Cross-normalisation
           case_control_norm_dfs[[1]][[1]], # scenario 1, cross-normalised (JYNR210602)
           case_control_norm_dfs[[5]][[1]], # scenario 2, cross-normalised (JYNR210602)
           limma_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Limma (1-4)
           case_control_norm_dfs[[3]][[1]], # scenario 1, limma (JYNR210602)
           case_control_norm_dfs[[7]][[1]], # scenario 2, limma (JYNR210602)
           cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 3,]), # Calibrator-normalised (indi SFs)
      list(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
           crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Cross-normalisation
           case_control_norm_dfs[[2]][[1]], # scenario 1, cross-normalised (JYNR220601)
           case_control_norm_dfs[[6]][[1]], # scenario 2, cross-normalised (JYNR220601)
           limma_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Limma (1-4)
           case_control_norm_dfs[[4]][[1]], # scenario 1, limma (JYNR220601)
           case_control_norm_dfs[[8]][[1]], # scenario 2, limma (JYNR220601)
           log_JYNR220601[meta_mouse_JYNR220601$Group == 2,])), #  Calibrator-normalised (indi SFs)
    placebos = list(
      list(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
           crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Cross-normalisation
           case_control_norm_dfs[[1]][[2]], # scenario 1, cross-normalised (JYNR210602)
           case_control_norm_dfs[[5]][[2]], # scenario 2, cross-normalised (JYNR210602)
           limma_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Limma (1-4)
           case_control_norm_dfs[[3]][[2]], # scenario 1, limma (JYNR210602)
           case_control_norm_dfs[[7]][[2]], # scenario 2, limma (JYNR210602)
           log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]), #  Calibrator-normalised (indi SFs)
      list(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # No normalisation (1-4)
           crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Cross-normalisation
           case_control_norm_dfs[[2]][[2]], # scenario 1, cross-normalised (JYNR220601)
           case_control_norm_dfs[[6]][[2]], # scenario 2, cross-normalised (JYNR220601)
           limma_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Limma (1-4)
           case_control_norm_dfs[[4]][[2]], # scenario 2, limma (JYNR220601)
           case_control_norm_dfs[[8]][[2]], # scenario 2, limma (JYNR220601)
           cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 1,])), # Calibrator-normalised (indi SFs)
    norm_methods = c("Baseline", "Cross", "Cross, real 1", "Cross, real 2",
                     "limma", "limma, real 1", "limma, real 2", "Cal"),
    short_y_lab = TRUE, fixed_ref = TRUE, #just_numbers = TRUE,
    x_angle = 30, x_hjust = 1, x_vjust = 1, meta_aligned = TRUE)
  
  # High quality filter
  pre_plot_h <- precision_plot(
    cases = list(
      list(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # No normalisation (1-4)
           crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Cross-normalisation
           case_control_norm_dfs[[1]][[1]], # scenario 1, cross-normalised (JYNR210602)
           case_control_norm_dfs[[5]][[1]], # scenario 2, cross-normalised (JYNR210602)
           limma_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Limma (1-4)
           case_control_norm_dfs[[3]][[1]], # scenario 1, limma (JYNR210602)
           case_control_norm_dfs[[7]][[1]], # scenario 2, limma (JYNR210602)
           cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 3,]), # Calibrator-normalised (indi SFs)
      list(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
           crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Cross-normalisation
           case_control_norm_dfs[[2]][[1]], # scenario 1, cross-normalised (JYNR220601)
           case_control_norm_dfs[[6]][[1]], # scenario 2, cross-normalised (JYNR220601)
           limma_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Limma (1-4)
           case_control_norm_dfs[[4]][[1]], # scenario 1, limma (JYNR220601)
           case_control_norm_dfs[[8]][[1]], # scenario 2, limma (JYNR220601)
           log_JYNR220601[meta_mouse_JYNR220601$Group == 2,])), #  Calibrator-normalised (indi SFs)
    placebos = list(
      list(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
           crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Cross-normalisation
           case_control_norm_dfs[[1]][[2]], # scenario 1, cross-normalised (JYNR210602)
           case_control_norm_dfs[[5]][[2]], # scenario 2, cross-normalised (JYNR210602)
           limma_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Limma (1-4)
           case_control_norm_dfs[[3]][[2]], # scenario 1, limma (JYNR210602)
           case_control_norm_dfs[[7]][[2]], # scenario 2, limma (JYNR210602)
           log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]), #  Calibrator-normalised (indi SFs)
      list(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # No normalisation (1-4)
           crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Cross-normalisation
           case_control_norm_dfs[[2]][[2]], # scenario 1, cross-normalised (JYNR220601)
           case_control_norm_dfs[[6]][[2]], # scenario 2, cross-normalised (JYNR220601)
           limma_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Limma (1-4)
           case_control_norm_dfs[[4]][[2]], # scenario 2, limma (JYNR220601)
           case_control_norm_dfs[[8]][[2]], # scenario 2, limma (JYNR220601)
           cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 1,])), # Calibrator-normalised (indi SFs)
    norm_methods = c("Baseline", "Cross", "Cross, real 1", "Cross, real 2",
                     "limma", "limma, real 1", "limma, real 2", "Cal"),
    short_y_lab = TRUE, fixed_ref = TRUE, #just_numbers = TRUE,
    x_angle = 30, x_hjust = 1, x_vjust = 1, meta_aligned = TRUE, apt_qual = "high")
  
  # Medium quality filter
  pre_plot_m <- precision_plot(
    cases = list(
      list(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # No normalisation (1-4)
           crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Cross-normalisation
           case_control_norm_dfs[[1]][[1]], # scenario 1, cross-normalised (JYNR210602)
           case_control_norm_dfs[[5]][[1]], # scenario 2, cross-normalised (JYNR210602)
           limma_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Limma (1-4)
           case_control_norm_dfs[[3]][[1]], # scenario 1, limma (JYNR210602)
           case_control_norm_dfs[[7]][[1]], # scenario 2, limma (JYNR210602)
           cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 3,]), # Calibrator-normalised (indi SFs)
      list(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
           crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Cross-normalisation
           case_control_norm_dfs[[2]][[1]], # scenario 1, cross-normalised (JYNR220601)
           case_control_norm_dfs[[6]][[1]], # scenario 2, cross-normalised (JYNR220601)
           limma_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Limma (1-4)
           case_control_norm_dfs[[4]][[1]], # scenario 1, limma (JYNR220601)
           case_control_norm_dfs[[8]][[1]], # scenario 2, limma (JYNR220601)
           log_JYNR220601[meta_mouse_JYNR220601$Group == 2,])), #  Calibrator-normalised (indi SFs)
    placebos = list(
      list(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
           crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Cross-normalisation
           case_control_norm_dfs[[1]][[2]], # scenario 1, cross-normalised (JYNR210602)
           case_control_norm_dfs[[5]][[2]], # scenario 2, cross-normalised (JYNR210602)
           limma_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Limma (1-4)
           case_control_norm_dfs[[3]][[2]], # scenario 1, limma (JYNR210602)
           case_control_norm_dfs[[7]][[2]], # scenario 2, limma (JYNR210602)
           log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]), #  Calibrator-normalised (indi SFs)
      list(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # No normalisation (1-4)
           crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Cross-normalisation
           case_control_norm_dfs[[2]][[2]], # scenario 1, cross-normalised (JYNR220601)
           case_control_norm_dfs[[6]][[2]], # scenario 2, cross-normalised (JYNR220601)
           limma_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Limma (1-4)
           case_control_norm_dfs[[4]][[2]], # scenario 2, limma (JYNR220601)
           case_control_norm_dfs[[8]][[2]], # scenario 2, limma (JYNR220601)
           cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 1,])), # Calibrator-normalised (indi SFs)
    norm_methods = c("Baseline", "Cross", "Cross, real 1", "Cross, real 2",
                     "limma", "limma, real 1", "limma, real 2", "Cal"),
    short_y_lab = TRUE, fixed_ref = TRUE,
    x_angle = 30, x_hjust = 1, x_vjust = 1, meta_aligned = TRUE, apt_qual = "med")
  
  # Low quality filter
  pre_plot_l <- precision_plot(
    cases = list(
      list(log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # No normalisation (1-4)
           crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Cross-normalisation
           case_control_norm_dfs[[1]][[1]], # scenario 1, cross-normalised (JYNR210602)
           case_control_norm_dfs[[5]][[1]], # scenario 2, cross-normalised (JYNR210602)
           limma_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Limma (1-4)
           case_control_norm_dfs[[3]][[1]], # scenario 1, limma (JYNR210602)
           case_control_norm_dfs[[7]][[1]], # scenario 2, limma (JYNR210602)
           cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 3,]), # Calibrator-normalised (indi SFs)
      list(log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
           crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Cross-normalisation
           case_control_norm_dfs[[2]][[1]], # scenario 1, cross-normalised (JYNR220601)
           case_control_norm_dfs[[6]][[1]], # scenario 2, cross-normalised (JYNR220601)
           limma_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Limma (1-4)
           case_control_norm_dfs[[4]][[1]], # scenario 1, limma (JYNR220601)
           case_control_norm_dfs[[8]][[1]], # scenario 2, limma (JYNR220601)
           log_JYNR220601[meta_mouse_JYNR220601$Group == 2,])), #  Calibrator-normalised (indi SFs)
    placebos = list(
      list(log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
           crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Cross-normalisation
           case_control_norm_dfs[[1]][[2]], # scenario 1, cross-normalised (JYNR210602)
           case_control_norm_dfs[[5]][[2]], # scenario 2, cross-normalised (JYNR210602)
           limma_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Limma (1-4)
           case_control_norm_dfs[[3]][[2]], # scenario 1, limma (JYNR210602)
           case_control_norm_dfs[[7]][[2]], # scenario 2, limma (JYNR210602)
           log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]), #  Calibrator-normalised (indi SFs)
      list(log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # No normalisation (1-4)
           crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Cross-normalisation
           case_control_norm_dfs[[2]][[2]], # scenario 1, cross-normalised (JYNR220601)
           case_control_norm_dfs[[6]][[2]], # scenario 2, cross-normalised (JYNR220601)
           limma_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Limma (1-4)
           case_control_norm_dfs[[4]][[2]], # scenario 2, limma (JYNR220601)
           case_control_norm_dfs[[8]][[2]], # scenario 2, limma (JYNR220601)
           cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 1,])), # Calibrator-normalised (indi SFs)
    norm_methods = c("Baseline", "Cross", "Cross, real 1", "Cross, real 2",
                     "limma", "limma, real 1", "limma, real 2", "Cal"),
    short_y_lab = TRUE, fixed_ref = TRUE,
    x_angle = 30, x_hjust = 1, x_vjust = 1, meta_aligned = TRUE, apt_qual = "low")
  
  ## Putting it all together
  sup_fig_13 <- ggarrange(pre_plot_b,NULL,pre_plot_h,NULL,pre_plot_m,NULL,pre_plot_l,
                          nrow = 7, heights = c(1,0.05,1,0.05,1,0.05,1),
                          labels = c("A", "", "B", "", "C", "", "D"))
  
  
  ## Save it
  fact <- 1.5
  ggsave("figures/supfig13.png", plot = sup_fig_13, width = fact*(210-20), height = fact*(297-30-25),
         units = "mm", dpi = 185, bg = "white")
}


################################################################################
## Supplementary figure 14: FC/alpha optimisation (baseline)
for (generate_figure in "sup_fig_14"){
  
  # Get references for both alpha and FC tests
  refs <- list(get_signif_targets(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                                  group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = 0.3),
               get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                                  group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = 0.3))
  
  ## Test different alpha levels
  for (test_alpha in "yes"){
    alphas <- c(0.05, 10^-2, 10^-3, 10^-4, 10^-5)
    fc_cutoff <- 0.3
    combinations <- length(alphas) * 2
    alpha_data <- data.frame(dataset = c(rep("SomaMouse1", combinations/2), rep("SomaMouse2", combinations/2)),
                             alpha = rep(NA, combinations),
                             precision = rep(NA, combinations),
                             positives = rep(NA, combinations))
    cases <- c("JYNR210602", "JYNR220601"); controls <- rev(cases)
    df_names <- c("SomaMouse1", "SomaMouse2")
    metas <- list(meta_mouse_JYNR210602, meta_mouse_JYNR220601)
    group1 <- c(3,2); contr <- 1
    
    i <- 0
    for (j in 1:2){
      dataset <- df_names[j]; print(dataset)
      ref <- refs[[j]]
      case <- group1[j]
      meta1 <- metas[[j]]; meta2 <- rev(metas)[[j]]
      soma1 <- get(paste0("log_", cases[j])); soma2 <- get(paste0("log_", controls[j]))
      
      for(alpha in alphas){
        i <- i + 1
        print(i)
        
        positive <- get_signif_targets(soma1 = soma1, soma2 = soma2,
                                       meta1 = meta1, meta2 = meta2,
                                       group1 = case, group2 = contr,
                                       fc_cutoff = fc_cutoff, alpha = alpha,
                                       sig_list_only = TRUE)
        TP <- sum(ref %in% positive)
        FP <- sum(!positive %in% ref)
        precision <- TP/(TP+FP)
        
        alpha_data[i,] <- c(dataset,
                            alpha,
                            precision,
                            length(positive))
      }
    }
  }
  
  ## Test different FC thresholds
  for (test_fc in "yes"){
    alpha <- 0.05
    fc_cutoffs <- c(0.3,0.6,0.9,1.2,1.5)
    combinations <- length(fc_cutoffs) * 2
    fc_data <- data.frame(dataset = c(rep("SomaMouse1", combinations/2), rep("SomaMouse2", combinations/2)),
                          FC_cutoff = rep(NA, combinations),
                          precision = rep(NA, combinations),
                          positives = rep(NA, combinations))
    cases <- c("JYNR210602", "JYNR220601"); controls <- rev(cases)
    df_names <- c("SomaMouse1", "SomaMouse2")
    metas <- list(meta_mouse_JYNR210602, meta_mouse_JYNR220601)
    group1 <- c(3,2); contr <- 1
    
    i <- 0
    for (j in 1:2){
      dataset <- df_names[j]; print(dataset)
      ref <- refs[[j]]
      case <- group1[j]
      meta1 <- metas[[j]]; meta2 <- rev(metas)[[j]]
      soma1 <- get(paste0("log_", cases[j])); soma2 <- get(paste0("log_", controls[j]))
      
      for(fc_cutoff in fc_cutoffs){
        i <- i + 1
        print(i)
        
        positive <- get_signif_targets(soma1 = soma1, soma2 = soma2,
                                       meta1 = meta1, meta2 = meta2,
                                       group1 = case, group2 = contr,
                                       fc_cutoff = fc_cutoff, alpha = alpha,
                                       sig_list_only = TRUE)
        TP <- sum(ref %in% positive)
        FP <- sum(!positive %in% ref)
        precision <- TP/(TP+FP)
        
        fc_data[i,] <- c(dataset,
                         fc_cutoff,
                         precision,
                         length(positive))
      }
    }
  }
  
  ## Pre-process
  # Alpha
  alpha_data[,2:4] <- lapply(2:4, function(x) as.numeric(alpha_data[,x]))
  l_a <- dim(alpha_data)[1]/2
  alpha_data$med_precision <- sapply(1:l_a,
                                     function(i) median(c(alpha_data$precision[i], alpha_data$precision[i+l_a])))
  alpha_data$med_pos <- sapply(1:l_a,
                                     function(i) median(c(alpha_data$positives[i], alpha_data$positives[i+l_a])))
  
  # FC
  fc_data[,2:4] <- lapply(2:4, function(x) as.numeric(fc_data[,x]))
  l_f <- dim(fc_data)[1]/2
  fc_data$med_precision <- sapply(1:l_f,
                                  function(i) median(c(fc_data$precision[i], fc_data$precision[i+l_a])))
  fc_data$med_pos <- sapply(1:l_f,
                                  function(i) median(c(fc_data$positives[i], fc_data$positives[i+l_a])))
  
  ## Plot it
  # Alpha, precision
  sup_fig_14a_1 <- alpha_data %>% ggplot(mapping = aes(x = alpha,
                                                       y = precision,
                                                       colour = dataset)) +
    ylim(c(0,1)) +
    geom_line() +
    geom_line(aes(x = alpha, y = med_precision), linetype = "dashed",
              colour = "#78ED56") +
    scale_x_reverse() +
    xlab("Significance threshold (alpha)") +
    ylab("Precision") +
    labs(colour = "Dataset") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.direction = "horizontal")   +
    annotation_compass("log2(FC) threshold = 0.3", "SW") +
    scale_colour_manual(values = c("#FFA355", "#55B1FF"))
  
  # Alpha, DE targets
  x1_start <- min(alpha_data$alpha); x1_end <- max(alpha_data$alpha)
  ymax1 <- max(alpha_data$positives)
  sup_fig_14a_2 <- alpha_data %>% ggplot(mapping = aes(x = alpha,
                                                       y = positives,
                                                       color = dataset)) +
    geom_line() +
    geom_line(aes(x = alpha, y = med_pos), linetype = "dashed",
              colour = "#78ED56") +
    xlab("Significance threshold (alpha)") +
    ylab("DE targets after swapping") +
    labs(colour = "Dataset") +
    theme_classic() +
    ylim(0,ymax1) +
    theme(plot.title = element_text(hjust = 0.5))  +
    annotation_compass("log2(FC) threshold = 0.3", "NE") +
    geom_segment(x = -0, xend = -x1_end, y = 100, yend = 100,
                 linetype = "dashed", colour = "lightgrey") +
    geom_segment(x = 0, xend = -x1_end, y = 50, yend = 50,
                 linetype = "dashed", colour = "darkgrey") +
    geom_segment(x = 0, xend = -x1_end, y = 0, yend = 0,
                 linetype = "dashed", colour = "black") +
    scale_x_reverse() +
    scale_colour_manual(values = c("#FFA355", "#55B1FF"))
  
  sup_fig_14a <- ggarrange(sup_fig_14a_1, sup_fig_14a_2,
                           legend = "none")
  
  
  # FC, precision
  sup_fig_14b_1 <- fc_data %>% ggplot(mapping = aes(x = FC_cutoff,
                                                    y = precision,
                                                    colour = dataset)) +
    ylim(c(0,1)) +
    geom_line() +
    geom_line(aes(x = FC_cutoff, y = med_precision), linetype = "dashed",
              colour = "#78ED56") +
    xlab("Log2(fold change) threshold") +
    ylab("Precision") +
    labs(colour = "Dataset") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))   +
    annotation_compass("alpha = 0.05", "SW") +
    scale_colour_manual(values = c("#FFA355", "#55B1FF"))
  
  # Alpha, DE targets
  x1_start <- min(fc_data$FC_cutoff); x1_end <- max(fc_data$FC_cutoff)
  ymax1 <- max(fc_data$positives)
  sup_fig_14b_2 <- fc_data %>% ggplot(mapping = aes(x = FC_cutoff,
                                                    y = positives,
                                                    color = dataset)) +
    geom_line() +
    geom_line(aes(x = FC_cutoff, y = med_pos), linetype = "dashed",
              colour = "#78ED56") +
    xlab("Log2(fold change) threshold") +
    ylab("DE targets after swapping") +
    labs(colour = "Dataset") +
    theme_classic() +
    ylim(0,ymax1) +
    theme(plot.title = element_text(hjust = 0.5))  +
    annotation_compass("alpha = 0.05", "NE") +
    geom_segment(x = x1_start, xend = x1_end, y = 100, yend = 100,
                 linetype = "dashed", colour = "lightgrey") +
    geom_segment(x = x1_start, xend = x1_end, y = 50, yend = 50,
                 linetype = "dashed", colour = "darkgrey") +
    geom_segment(x = x1_start, xend = x1_end, y = 0, yend = 0,
                 linetype = "dashed", colour = "black") +
    scale_colour_manual(values = c("#FFA355", "#55B1FF"))
  
  sup_fig_14b <- ggarrange(sup_fig_14b_1, sup_fig_14b_2,
                           legend = "none")
  alpha_data %>% ggplot(mapping = aes(x = alpha,
                                      y = positives,
                                      color = dataset)) + geom_point()
  legend <- get_legend(sup_fig_14a_1)
  sup_fig_14 <- ggarrange(sup_fig_14a, sup_fig_14b,
            nrow = 2, labels = c("A", "B"),
            common.legend = TRUE,
            legend.grob = legend,
            legend = "top")
  
  fact <- 1.5
  ggsave("figures/supfig14.png", plot = sup_fig_14, width = fact*(210-20), height = fact*(297-30-180),
         units = "mm", dpi = 185, bg = "white")
}


################################################################################
## Supplementary figure 15: FC/alpha optimisation (high_qual aptamers)
for (generate_figure in "sup_fig_15"){
  
  # Set high quality SOMAmer filter
  seqids <- SomaScanAnnotation$SeqId[SomaScanAnnotation$mouse_quality %in% "High" & SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)]
  
  # Get references for both alpha and FC tests
  refs <- list(get_signif_targets(soma1 = log_JYNR210602[,colnames(log_JYNR210602) %in% seqids],
                                  meta1 = meta_mouse_JYNR210602,
                                  group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = 0.3),
               get_signif_targets(soma1 = log_JYNR220601[,colnames(log_JYNR220601) %in% seqids],
                                  meta1 = meta_mouse_JYNR220601,
                                  group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = 0.3))
  
  ## Test different alpha levels
  for (test_alpha in "yes"){
    alphas <- c(0.05, 10^-2, 10^-3, 10^-4, 10^-5)
    fc_cutoff <- 0.3
    combinations <- length(alphas) * 2
    alpha_data <- data.frame(dataset = c(rep("SomaMouse1", combinations/2), rep("SomaMouse2", combinations/2)),
                             alpha = rep(NA, combinations),
                             precision = rep(NA, combinations),
                             positives = rep(NA, combinations))
    cases <- c("JYNR210602", "JYNR220601"); controls <- rev(cases)
    df_names <- c("SomaMouse1", "SomaMouse2")
    metas <- list(meta_mouse_JYNR210602, meta_mouse_JYNR220601)
    group1 <- c(3,2); contr <- 1
    
    i <- 0
    for (j in 1:2){
      dataset <- df_names[j]; print(dataset)
      ref <- refs[[j]]
      case <- group1[j]
      meta1 <- metas[[j]]; meta2 <- rev(metas)[[j]]
      soma1 <- get(paste0("log_", cases[j])); soma2 <- get(paste0("log_", controls[j]))
      
      for(alpha in alphas){
        i <- i + 1
        print(i)
        
        positive <- get_signif_targets(soma1 = soma1[,colnames(soma1) %in% seqids],
                                       soma2 = soma2[,colnames(soma2) %in% seqids],
                                       meta1 = meta1, meta2 = meta2,
                                       group1 = case, group2 = contr,
                                       fc_cutoff = fc_cutoff, alpha = alpha,
                                       sig_list_only = TRUE)
        TP <- sum(ref %in% positive)
        FP <- sum(!positive %in% ref)
        precision <- TP/(TP+FP)
        
        alpha_data[i,] <- c(dataset,
                            alpha,
                            precision,
                            length(positive))
      }
    }
  }
  
  ## Test different FC thresholds
  for (test_fc in "yes"){
    alpha <- 0.05
    fc_cutoffs <- c(0.3,0.6,0.9,1.2,1.5)
    combinations <- length(fc_cutoffs) * 2
    fc_data <- data.frame(dataset = c(rep("SomaMouse1", combinations/2), rep("SomaMouse2", combinations/2)),
                          FC_cutoff = rep(NA, combinations),
                          precision = rep(NA, combinations),
                          positives = rep(NA, combinations))
    cases <- c("JYNR210602", "JYNR220601"); controls <- rev(cases)
    df_names <- c("SomaMouse1", "SomaMouse2")
    metas <- list(meta_mouse_JYNR210602, meta_mouse_JYNR220601)
    group1 <- c(3,2); contr <- 1
    
    i <- 0
    for (j in 1:2){
      dataset <- df_names[j]; print(dataset)
      ref <- refs[[j]]
      case <- group1[j]
      meta1 <- metas[[j]]; meta2 <- rev(metas)[[j]]
      soma1 <- get(paste0("log_", cases[j])); soma2 <- get(paste0("log_", controls[j]))
      
      for(fc_cutoff in fc_cutoffs){
        i <- i + 1
        print(i)
        
        positive <- get_signif_targets(soma1 = soma1[,colnames(soma1) %in% seqids],
                                       soma2 = soma2[,colnames(soma2) %in% seqids],
                                       meta1 = meta1, meta2 = meta2,
                                       group1 = case, group2 = contr,
                                       fc_cutoff = fc_cutoff, alpha = alpha,
                                       sig_list_only = TRUE)
        TP <- sum(ref %in% positive)
        FP <- sum(!positive %in% ref)
        precision <- TP/(TP+FP)
        
        fc_data[i,] <- c(dataset,
                         fc_cutoff,
                         precision,
                         length(positive))
      }
    }
  }
  
  ## Pre-process
  # Alpha
  alpha_data[,2:4] <- lapply(2:4, function(x) as.numeric(alpha_data[,x]))
  l_a <- dim(alpha_data)[1]/2
  alpha_data$med_precision <- sapply(1:l_a,
                                     function(i) median(c(alpha_data$precision[i], alpha_data$precision[i+l_a])))
  alpha_data$med_pos <- sapply(1:l_a,
                               function(i) median(c(alpha_data$positives[i], alpha_data$positives[i+l_a])))
  
  # FC
  fc_data[,2:4] <- lapply(2:4, function(x) as.numeric(fc_data[,x]))
  l_f <- dim(fc_data)[1]/2
  fc_data$med_precision <- sapply(1:l_f,
                                  function(i) median(c(fc_data$precision[i], fc_data$precision[i+l_a])))
  fc_data$med_pos <- sapply(1:l_f,
                            function(i) median(c(fc_data$positives[i], fc_data$positives[i+l_a])))
  
  ## Plot it
  # Alpha, precision
  sup_fig_15a_1 <- alpha_data %>% ggplot(mapping = aes(x = alpha,
                                                       y = precision,
                                                       colour = dataset)) +
    ylim(c(0,1)) +
    geom_line() +
    geom_line(aes(x = alpha, y = med_precision), linetype = "dashed",
              colour = "#78ED56") +
    scale_x_reverse() +
    xlab("Significance threshold (alpha)") +
    ylab("Precision") +
    labs(colour = "Dataset") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.direction = "horizontal")   +
    annotation_compass("log2(FC) threshold = 0.3", "SW") +
    scale_colour_manual(values = c("#FFA355", "#55B1FF"))
  
  # Alpha, DE targets
  x1_start <- min(alpha_data$alpha); x1_end <- max(alpha_data$alpha)
  ymax1 <- max(alpha_data$positives)
  sup_fig_15a_2 <- alpha_data %>% ggplot(mapping = aes(x = alpha,
                                                       y = positives,
                                                       color = dataset)) +
    geom_line() +
    geom_line(aes(x = alpha, y = med_pos), linetype = "dashed",
              colour = "#78ED56") +
    xlab("Significance threshold (alpha)") +
    ylab("DE targets after swapping") +
    labs(colour = "Dataset") +
    theme_classic() +
    ylim(0,ymax1) +
    theme(plot.title = element_text(hjust = 0.5))  +
    annotation_compass("log2(FC) threshold = 0.3", "NE") +
    geom_segment(x = -0, xend = -x1_end, y = 100, yend = 100,
                 linetype = "dashed", colour = "lightgrey") +
    geom_segment(x = 0, xend = -x1_end, y = 50, yend = 50,
                 linetype = "dashed", colour = "darkgrey") +
    geom_segment(x = 0, xend = -x1_end, y = 0, yend = 0,
                 linetype = "dashed", colour = "black") +
    scale_x_reverse() +
    scale_colour_manual(values = c("#FFA355", "#55B1FF"))
  
  sup_fig_15a <- ggarrange(sup_fig_15a_1, sup_fig_15a_2,
                           legend = "none")
  
  
  # FC, precision
  sup_fig_15b_1 <- fc_data %>% ggplot(mapping = aes(x = FC_cutoff,
                                                    y = precision,
                                                    colour = dataset)) +
    ylim(c(0,1)) +
    geom_line() +
    geom_line(aes(x = FC_cutoff, y = med_precision), linetype = "dashed",
              colour = "#78ED56") +
    xlab("Log2(fold change) threshold") +
    ylab("Precision") +
    labs(colour = "Dataset") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))   +
    annotation_compass("alpha = 0.05", "SW") +
    scale_colour_manual(values = c("#FFA355", "#55B1FF"))
  
  # FC, DE targets
  x1_start <- min(fc_data$FC_cutoff); x1_end <- max(fc_data$FC_cutoff)
  ymax1 <- max(fc_data$positives)
  sup_fig_15b_2 <- fc_data %>% ggplot(mapping = aes(x = FC_cutoff,
                                                    y = positives,
                                                    color = dataset)) +
    geom_line() +
    geom_line(aes(x = FC_cutoff, y = med_pos), linetype = "dashed",
              colour = "#78ED56") +
    xlab("Log2(fold change) threshold") +
    ylab("DE targets after swapping") +
    labs(colour = "Dataset") +
    theme_classic() +
    ylim(0,ymax1) +
    theme(plot.title = element_text(hjust = 0.5))  +
    annotation_compass("alpha = 0.05", "NE") +
    geom_segment(x = x1_start, xend = x1_end, y = 100, yend = 100,
                 linetype = "dashed", colour = "lightgrey") +
    geom_segment(x = x1_start, xend = x1_end, y = 50, yend = 50,
                 linetype = "dashed", colour = "darkgrey") +
    geom_segment(x = x1_start, xend = x1_end, y = 0, yend = 0,
                 linetype = "dashed", colour = "black") +
    scale_colour_manual(values = c("#FFA355", "#55B1FF"))
  
  sup_fig_15b <- ggarrange(sup_fig_15b_1, sup_fig_15b_2,
                           legend = "none")
  alpha_data %>% ggplot(mapping = aes(x = alpha,
                                      y = positives,
                                      color = dataset)) + geom_point()
  legend <- get_legend(sup_fig_15a_1)
  sup_fig_15 <- ggarrange(sup_fig_15a, sup_fig_15b,
                          nrow = 2, labels = c("A", "B"),
                          common.legend = TRUE,
                          legend.grob = legend,
                          legend = "top")
  
  fact <- 1.5
  ggsave("figures/supfig15.png", plot = sup_fig_15, width = fact*(210-20), height = fact*(297-30-180),
         units = "mm", dpi = 185, bg = "white")
}

# high_fc_data <- fc_data; high_alpha_data <- alpha_data








################################################################################
## Supplementary figure 16: Top 100 intra-study overlap (boostrapping)
for (generate_figure in "sup_fig_16"){
  # Get data
  sup_fig_16_data <- bootstrap_data_supfig16
  
  # Compute plot locations
  y1 <- median(sup_fig_16_data$overlaps[sup_fig_16_data$dataset == "SomaMouse1"])
  y2 <- median(sup_fig_16_data$overlaps[sup_fig_16_data$dataset == "SomaMouse2"])
  x1_start <- 0.5;
  x1_end <- dim(sup_fig_16_data)[1]/2 + 0.5;
  x2_start <- 0.5 + dim(sup_fig_16_data)[1]/2;
  x2_end <- dim(sup_fig_16_data)[1] + 0.5;
  
  # Create the plot
  sup_fig_16 <- sup_fig_16_data %>%
    ggplot(mapping = aes(x = iteration,
                         y = overlaps,
                         fill = dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.direction = "vertical") +
    annotate("text", label = paste0("Median = ", y1),
             x = x1_start, y = 100, hjust = 0) +
    annotate("text", label = paste0("Median = ", y2),
             x = x2_end, y = 100, hjust = 1) +
    scale_fill_manual(values = c("#FFA355", "#55B1FF")) +
    geom_segment(x = x1_start, xend = x1_end, y = y1, yend = y1,
                 linetype = "dashed", color = "#FFA355") +
    geom_segment(x = x2_start, xend = x2_end, y = y2, yend = y2,
                 linetype = "dashed", color = "#55B1FF") +
    xlab("Bootstrap iteration") +
    ylab("Overlaps between top100 genes") +
    labs(fill = "Dataset") +
    ggtitle("Top 100 SOMAmers by FC in original dataset versus boostrapped dataset")
  
  ggsave("figures/supfig16.png", plot = sup_fig_16, width = fact*(210-20), height = fact*(297-30-180),
         units = "mm", dpi = 185, bg = "white")
}
