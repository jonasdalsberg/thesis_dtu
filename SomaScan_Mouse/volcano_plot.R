## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel); library(grid)
setwd("~/JDJG/clean")


## Load data
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")
source("scripts/SomaScan_Mouse/volcano_function.R")

################################################################################
# Plots
################################################################################


# Run this chunk to get all standard mouse plots
# (DIO-NASH Vehicle vs. Chow Vehicle in both studies + control:control plots)
# (for all cross-normalization options)
for (get_all_mouse_plots in "go"){
  # Plot and save (designed for plotting the best contrasts in both studies)
  legend <- get_legend(volcano_plot(legend.dir = "horizontal"))

  for (norm in c(FALSE, "soma_cross", "limma")){
  # for (norm in c(FALSE)){
    for (oneclickboom in "Letsgo"){
      # Parameters
      # norm <- "MZCG_cross"
      # norm <- "limma"
      # norm <- FALSE
      alpha <- 0.05
      # alpha <- 0.1
      mtc <- "none"
      mtc <- "bh"
      
      if (norm == "MZCG_cross"){
        soma1 <- "norm_JYNR210602"; soma2 <- "norm_JYNR220601"
        title <- paste0(title, ", cross-normalised (MZCG)")
        print("MZCG cross-norm")
      } else if (norm == "soma_cross"){
        soma1 <- "crossnorm_JYNR210602"; soma2 <- "crossnorm_JYNR220601"
        title <- paste0(title, ", cross-normalised (SomaLogic)")
        print ("SomaLogic cross-norm")
      } else if (norm == "limma"){
        soma1 <- "limma_JYNR210602"; soma2 <- "limma_JYNR220601"
        title <- paste0(title, ", batch effects corrected (limma)")
        print ("limma")
      } else {
        soma1 <- "log_JYNR210602"; soma2 <- "log_JYNR220601"
        print("log")
      }
      soma1 <- get(paste0(soma1)); soma2 <- get(paste0(soma2))
      
      # Plot
      myplots <- list(volcano_plot(soma1 = soma1, alpha = alpha, mtc = mtc, title = "SomaMouse1"),
                      volcano_plot(soma1 = soma2, meta1 = meta_mouse_JYNR220601, title = "SomaMouse2",
                                   alpha = alpha, mtc = mtc,
                                   control = 1, case = 2))
      
      g <- ggarrange(plotlist = myplots,
                     ncol = length(myplots),
                     legend.grob = legend,
                     hjust = -2)
      title <- "DIO-NASH Vehicle vs. Chow Vehicle"
      filename <- "figures/JYNR210602_JYNR220601_volcano_sep"
      if (norm == "MZCG_cross"){
        title <- paste0(title, ", cross-normalised (MZCG)")
        filename <- paste0(filename, "_mzcg_crossnorm")
      }
      if (norm == "soma_cross"){
        title <- paste0(title, ", cross-normalised (SomaLogic)")
        filename <- paste0(filename, "_soma_crossnorm")
      }
      if (norm == "limma"){
        title <- paste0(title, ", batch effects corrected (limma)")
        filename <- paste0(filename, "_limma")
      }
      if (!mtc == "none"){
        filename <- paste0(filename, "_", mtc)
      }
      if (alpha != 0.05){
        filename <- paste0(filename, "_alph_", alpha)
      }
      t
      filename <- paste0(filename, ".png")
      ggsave(filename, width = 4000, height = 1000, units = "px", dpi = 185, bg = "white")
    }
  }
  
  
  # Plot and save (designed for plotting control:control differences btw. mouse studies)
  for (norm in c(FALSE, "soma_cross", "limma")){
  # for (norm in c(FALSE)){
    for (oneclickboom in "Letsgo"){
      # Parameters
      # norm <- "MZCG_cross"
      # norm <- "limma"
      # norm <- FALSE
      alpha <- 0.05
      # alpha <- 0.1
      mtc <- "none"
      mtc <- "bh"
      
      if (norm == "MZCG_cross"){
        soma1 <- "norm_JYNR210602"; soma2 <- "norm_JYNR220601"
        title <- paste0(title, ", cross-normalised (MZCG)")
        print("MZCG cross-norm")
      } else if (norm == "soma_cross"){
        soma1 <- "crossnorm_JYNR210602"; soma2 <- "crossnorm_JYNR220601"
        title <- paste0(title, ", cross-normalised (SomaLogic)")
        print ("SomaLogic cross-norm")
      } else if (norm == "limma"){
        soma1 <- "limma_JYNR210602"; soma2 <- "limma_JYNR220601"
        title <- paste0(title, ", batch effects corrected (limma)")
        print ("limma")
      } else {
        soma1 <- "log_JYNR210602"; soma2 <- "log_JYNR220601"
        print("log")
      }
      soma1 <- get(paste0(soma1)); soma2 <- get(paste0(soma2))
      
      
      # Plot
      myplots <- list(volcano_plot(soma1 = soma1, soma2 = soma2, meta2 = meta_mouse_JYNR220601,
                                   group1 = 1, group2 = 1, title = "Chow Vehicle",
                                   alpha = alpha, mtc = mtc),
                      volcano_plot(soma1 = soma1, soma2 = soma2, meta2 = meta_mouse_JYNR220601,
                                   group1 = 3, group2 = 2, title = "DIO-NASH Vehicle",
                                   alpha = alpha, mtc = mtc))
      
      g <- ggarrange(plotlist = myplots,
                     ncol = length(myplots),
                     legend.grob = legend,
                     hjust = -2)
      title <- "SomaMouse1 vs. SomaMouse2"
      filename <- "figures/JYNR210602_JYNR220601_volcano_contr"
      if (norm == "MZCG_cross"){
        title <- paste0(title, ", cross-normalised (MZCG)")
        filename <- paste0(filename, "_mzcg_crossnorm")
      }
      if (norm == "soma_cross"){
        title <- paste0(title, ", cross-normalised (SomaLogic)")
        filename <- paste0(filename, "_soma_crossnorm")
      }
      if (norm == "limma"){
        title <- paste0(title, ", batch effects corrected (limma)")
        filename <- paste0(filename, "_limma")
      }
      if (!mtc == "none"){
        filename <- paste0(filename, "_", mtc)
      }
      if (alpha != 0.05){
        filename <- paste0(filename, "_alph_", alpha)
      }
      
      filename <- paste0(filename, ".png")
      ggsave(filename, width = 4000, height = 1000, units = "px", dpi = 185, bg = "white")
    }
  }
}


# Run this chunk to get placebo-swapped contrasts
# (for all cross-normalization options)
for (get_all_mouse_plots in "go"){
  legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
  # one_to_four <- FALSE # comment to run with the dfs lacking the last step (mednormSMP)
  for (norm in c(FALSE, "soma_cross", "limma")){
    for (oneclickboom in "Letsgo"){
      # Parameters
      # norm <- TRUE
      # norm <- "limma"
      # norm <- FALSE
      alpha <- 0.05
      # alpha <- 0.1
      mtc <- "none"
      mtc <- "bh"
      
      if (norm == "MZCG_cross"){
        soma1 <- "norm_JYNR210602"; soma2 <- "norm_JYNR220601"
        title <- paste0(title, ", cross-normalised (MZCG)")
        print("MZCG cross-norm")
      } else if (norm == "soma_cross"){
        soma1 <- "crossnorm_JYNR210602"; soma2 <- "crossnorm_JYNR220601"
        title <- paste0(title, ", cross-normalised (SomaLogic)")
        print ("SomaLogic cross-norm")
      } else if (norm == "limma"){
        soma1 <- "limma_JYNR210602"; soma2 <- "limma_JYNR220601"
        title <- paste0(title, ", batch effects corrected (limma)")
        print ("limma")
      } else {
        soma1 <- "log_JYNR210602"; soma2 <- "log_JYNR220601"
        print("log")
      }
      soma1 <- get(paste0(soma1)); soma2 <- get(paste0(soma2))
      
      
      # Plot
      myplots <- list(volcano_plot(soma1 = soma1, soma2 = soma2, meta2 = meta_mouse_JYNR220601,
                                   group1 = 3, group2 = 1, title = "SomaMouse1",
                                   alpha = alpha, mtc = mtc),
                      volcano_plot(soma1 = soma2, soma2 = soma1,
                                   meta1 = meta_mouse_JYNR220601, meta2 = meta_mouse_JYNR210602,
                                   group1 = 2, group2 = 1, title = "SomaMouse2",
                                   alpha = alpha, mtc = mtc))
      
      g <- ggarrange(plotlist = myplots,
                     ncol = length(myplots),
                     legend.grob = legend,
                     hjust = -2)
      title <- "DIO-NASH Vehicle vs. Chow Vehicle (swapped)"
      filename <- "figures/JYNR210602_JYNR220601_volcano_swapped"
      if (norm == "MZCG_cross"){
        title <- paste0(title, ", cross-normalised (MZCG)")
        filename <- paste0(filename, "_mzcg_crossnorm")
      }
      if (norm == "soma_cross"){
        title <- paste0(title, ", cross-normalised (SomaLogic)")
        filename <- paste0(filename, "_soma_crossnorm")
      }
      if (norm == "limma"){
        title <- paste0(title, ", batch effects corrected (limma)")
        filename <- paste0(filename, "_limma")
      }
      if (!mtc == "none"){
        filename <- paste0(filename, "_", mtc)
      }
      if (alpha != 0.05){
        filename <- paste0(filename, "_alph_", alpha)
      }
      annotate_figure(g, top = text_grob(title))
      filename <- paste0(filename, ".png")
      ggsave(filename, width = 4000, height = 1000, units = "px", dpi = 185, bg = "white")
    }
  }
}




# ################################################################################
# # This chunk is designed to plot the case_control_norm_dfs list of normalised dfs
# for (get_all_case_control_norm_plots in "go"){
#   # Plot and save (designed for plotting the best contrasts in both studies)
#   legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
#   n_plots <- length(case_control_norm_dfs)/2
#   norms <- c(TRUE, "limma", TRUE, "limma", FALSE, "limma", FALSE, "limma")
#   steps_incl <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
#   scenario <- c(1,1,2,2,3,3,4,4)
#   for (j in 1:n_plots){
#     print(j)
#     i <- j*2 - 1
#     # Parameters
#     # norm <- TRUE
#     # norm <- "limma"
#     # norm <- FALSE
#     alpha <- 0.05
#     # alpha <- 0.1
#     mtc <- "none"
#     mtc <- "bh"
#     
#     # Choose data
#     case1 <- case_control_norm_dfs[[i]][[1]]; placebo1 <- case_control_norm_dfs[[i+1]][[2]]
#     case2 <- case_control_norm_dfs[[i+1]][[1]]; placebo2 <- case_control_norm_dfs[[i]][[2]]
#     norm <- norms[j]; one_to_four <- steps_incl[j]
#     print("Data chosen.")
#     
#     # Plot intra-study contrasts
#     myplots <- list(volcano_plot(soma1 = case1, soma2 = placebo1,
#                                  alpha = alpha, mtc = mtc,
#                                  title = "JYNR210602", meta_aligned = TRUE),
#                     volcano_plot(soma1 = case2, soma2 = placebo2,
#                                  meta1 = meta_mouse_JYNR220601,
#                                  alpha = alpha, mtc = mtc,
#                                  title = "JYNR220601", meta_aligned = TRUE))
#     print("Intra-study plot created.")
#     g <- ggarrange(plotlist = myplots,
#                    ncol = length(myplots),
#                    common.legend = TRUE,
#                    legend.grob = legend,
#                    hjust = -2)
#     title <- "DIO-NASH Vehicle vs. Chow Vehicle"
#     if (norm == TRUE){
#       title <- paste0(title, ", cross-normalised")
#     }
#     if (norm == "limma"){
#       title <- paste0(title, ", batch effects corrected (limma)")
#     }
#     if (one_to_four == TRUE){
#       title <- paste0(title, ", 1-4 steps only")
#     }
#     title <- paste0(title, ", scenario ", scenario[j])
#     annotate_figure(g, top = text_grob(title))
#     print("Intra-study plot refined.")
#     
#     filename <- "figures/JYNR210602_JYNR220601_volcano_sep"
#     if (norm == TRUE){
#       filename <- paste0(filename, "_crossnorm")
#     }
#     if (norm == "limma"){
#       filename <- paste0(filename, "_limma")
#     }
#     if (one_to_four == TRUE){
#       filename <- paste0(filename, "_1_4")
#     }
#     if (!mtc == "none"){
#       filename <- paste0(filename, "_", mtc)
#     }
#     if (alpha != 0.05){
#       filename <- paste0(filename, "_alph_", alpha)
#     }
#     filename <- (paste0(filename, "_sce", scenario[j]))
#     filename <- paste0(filename, ".png")
#     ggsave(filename, width = 4000, height = 1200, units = "px", dpi = 200, bg = "white")
#     print("Intra-study plot saved.")
#     
#     
#     # Plot inter-study contrasts
#     myplots <- list(volcano_plot(soma1 = placebo1, soma2 = placebo2,
#                                  title = "Chow Vehicle", alpha = alpha,
#                                  mtc = mtc, meta_aligned = TRUE),
#                     volcano_plot(soma1 = case1, soma2 = case2,
#                                  title = "DIO-NASH Vehicle", alpha = alpha,
#                                  mtc = mtc, meta_aligned = TRUE))
#     print("Inter-study plot created.")
#     
#     g <- ggarrange(plotlist = myplots,
#                    ncol = length(myplots),
#                    common.legend = TRUE,
#                    legend.grob = legend,
#                    hjust = -2)
#     title <- "JYNR210602 vs. JYNR220601"
#     if (norm == TRUE){
#       title <- paste0(title, ", cross-normalised")
#     }
#     if (norm == "limma"){
#       title <- paste0(title, ", batch effects corrected (limma)")
#     }
#     if (one_to_four == TRUE){
#       title <- paste0(title, ", 1-4 steps only")
#     }
#     title <- paste0(title, ", scenario ", scenario[j])
#     annotate_figure(g, top = text_grob(title))
#     print("Inter-study plot refined.")
#     
#     filename <- "figures/JYNR210602_JYNR220601_volcano_contr"
#     if (norm == TRUE){
#       filename <- paste0(filename, "_crossnorm")
#     }
#     if (norm == "limma"){
#       filename <- paste0(filename, "_limma")
#     }
#     if (one_to_four == TRUE){
#       filename <- paste0(filename, "_1_4")
#     }
#     if (!mtc == "none"){
#       filename <- paste0(filename, "_", mtc)
#     }
#     if (alpha != 0.05){
#       filename <- paste0(filename, "_alph_", alpha)
#     }
#     filename <- (paste0(filename, "_sce", scenario[j]))
#     filename <- paste0(filename, ".png")
#     ggsave(filename, width = 4000, height = 1200, units = "px", dpi = 200, bg = "white")
#     print("Inter-study plot saved.")
#   }
# }



# ################################################################################
# # This chunk is designed to plot the case_control_norm_dfs list of normalised dfs
# # For PLACEBO SWAPPED plots
# for (get_all_case_control_norm_placebo_swapped_plots in "go"){
#   # Plot and save (designed for plotting the best contrasts in both studies)
#   legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
#   n_plots <- length(case_control_norm_dfs)/2 / 2
#   norms <- c(TRUE, "limma", TRUE, "limma", FALSE, "limma", FALSE, "limma")
#   steps_incl <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
#   scenario <- c(1,1,2,2,3,3,4,4)
#   for (j in 1:n_plots){
#     print(j)
#     i <- j*2 - 1
#     # Parameters
#     # norm <- TRUE
#     # norm <- "limma"
#     # norm <- FALSE
#     alpha <- 0.05
#     # alpha <- 0.1
#     mtc <- "none"
#     mtc <- "bh"
#     
#     # Choose data
#     case1 <- case_control_norm_dfs[[i]][[1]]; placebo1 <- case_control_norm_dfs[[i+1]][[2]]
#     case2 <- case_control_norm_dfs[[i+1]][[1]]; placebo2 <- case_control_norm_dfs[[i]][[2]]
#     norm <- norms[j]; one_to_four <- steps_incl[j]
#     print("Data chosen.")
#     
#     # Plot intra-study contrasts
#     myplots <- list(volcano_plot(soma1 = case1, soma2 = placebo2,
#                                  alpha = alpha, mtc = mtc,
#                                  title = "JYNR210602", meta_aligned = TRUE),
#                     volcano_plot(soma1 = case2, soma2 = placebo1,
#                                  alpha = alpha, mtc = mtc,
#                                  title = "JYNR220601", meta_aligned = TRUE))
#     print("Intra-study plot created (placebo swapped).")
#     
#     g <- ggarrange(plotlist = myplots,
#                    ncol = length(myplots),
#                    common.legend = TRUE,
#                    legend.grob = legend,
#                    hjust = -2)
#     title <- "DIO-NASH Vehicle vs. Chow Vehicle (swapped)"
#     filename <- "figures/JYNR210602_JYNR220601_volcano_swapped"
#     if (norm == TRUE){
#       title <- paste0(title, ", cross-normalised")
#       filename <- paste0(filename, "_crossnorm")
#     }
#     if (norm == "limma"){
#       title <- paste0(title, ", batch effects corrected (limma)")
#       filename <- paste0(filename, "_limma")
#     }
#     if (one_to_four == TRUE & (!norm == TRUE)){
#       title <- paste0(title, ", 1-4 steps only")
#       filename <- paste0(filename, "_1_4")
#     }
#     if (!mtc == "none"){
#       filename <- paste0(filename, "_", mtc)
#     }
#     if (alpha != 0.05){
#       filename <- paste0(filename, "_alph_", alpha)
#     }
#     annotate_figure(g, top = text_grob(title))
#     print("Intra-study plot refined (placebo swapped).")
#     filename <- paste0(filename, "_sce", scenario[j])
#     filename <- paste0(filename, ".png")
#     ggsave(filename, width = 4000, height = 1200, units = "px", dpi = 200, bg = "white")
#     print("Intra-study plot saved (placebo swapped).")
#   }
# }


# Plot and save
# legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
# for (oneclickboom in "Letsgo"){
#   # Choose groups and dataframe to be plotted
#   # 1 = JYNR210602, 2 = JYNR220601
#   df <- 1
#   norm <- TRUE
#   norm <- FALSE
#   alpha <- 0.05
#   # alpha <- 0.1
#   
#   
#   if (df == 1){
#     if (norm == TRUE){
#       soma_ <- norm_JYNR210602
#     } else {
#       soma_ <- log_JYNR210602
#     }
#     log_t_ <- TRUE
#     meta_ <- meta_mouse_JYNR210602
#     dataname <- "JYNR210602"
#     case <- 3
#     control <- 1
#   }
#   if (df == 2){
#     if (norm == TRUE){
#       soma_ <- norm_JYNR220601
#     } else {
#       soma_ <- log_JYNR220601
#     }
#     log_t_ <- TRUE
#     meta_ <- meta_mouse_JYNR220601
#     dataname <- "JYNR220601"
#     case <- 2
#     control <- 1
#   }
#   
#   
#   # Plot
#   myplots <- list(volcano_plot(soma1 = soma_, meta1 = meta_, title = FALSE, alpha = alpha,
#                                log_t = log_t_, control = control, case = case),
#                   volcano_plot(soma1 = soma_, meta1 = meta_, mtc = "BH", title = FALSE,
#                                alpha = alpha, log_t = log_t_, control = control, case = case),
#                   volcano_plot(soma1 = soma_, meta1 = meta_, mtc = "bonferroni", title = FALSE,
#                                alpha = alpha, log_t = log_t_, control = control, case = case))
#   
#   g <- ggarrange(plotlist = myplots,
#             ncol = length(myplots),
#             common.legend = TRUE,
#             legend.grob = legend,
#             hjust = -2)
#   annotate_figure(g, top = text_grob(paste0(dataname, ": Group ", case, " vs group ", control, " (control)")))
#   filename <- paste0("figures/", dataname, "_volcano")
#   filename <- paste0(filename, "_", case, "_v_", control)
#   if (norm == TRUE){
#     filename <- paste0(filename, "_crossnorm")
#   }
#   filename <- paste0(filename, "_mtc")
#   if (alpha != 0.05){
#     filename <- paste0(filename, "_alph_", alpha)
#   }
#   filename <- paste0(filename, ".png")
#   ggsave(filename, width = 3000, height = 1200, units = "px", dpi = 200, bg = "white")
# }


#################
# More customizable volcano plots (e.g. control:control)

# # Current (20/10/2023) version of control:control for mouse studies
# for (go in "go"){
#   norm = TRUE
#   norm = "limma"
#   norm = FALSE
#   mtc <- "none"
#   mtc <- "bh"
#   alpha <- 0.05
#   title = TRUE
#   legend.dir = "horizontal"
#   
#   if (norm == TRUE){
#     soma_1 <- norm_JYNR210602
#     soma_2 <- norm_JYNR220601
#   } else if (norm == "limma") {
#     soma_1 <- limma_JYNR210602
#     soma_2 <- limma_JYNR220601
#   } else {
#     soma_1 <- log_JYNR210602
#     soma_2 <- log_JYNR220601
#   }
#   soma <- soma_1
#   meta_1 <- meta_mouse_JYNR210602; meta_2 <- meta_mouse_JYNR220601
#   datanames <- c("JYNR210602", "JYNR220601")
#   myplots <- list()
#   grp_1_list <- c(1,3)
#   grp_2_list <- c(1,2)
#   
#   
#   
#   mtc <- stringr::str_to_lower(mtc)
#   mtc_title <- stringr::str_to_title(mtc)
#   if (mtc == "bh"){
#     mtc <- mtc_title <- "BH"
#   }
#   
#   for (i in 1:length(grp_1_list)){
#     grp_1 <- grp_1_list[i]
#     grp_2 <- grp_2_list[i]
#     # compute log2 fold change
#     fc <- sapply(1:ncol(soma_1), function(x) mean(na.omit(as.vector(unlist(soma_1[which(meta_1$Group==grp_1),x])))) - mean(na.omit(as.vector(unlist(soma_2[which(meta_2$Group==grp_2),x])))))
#     
#     # compute p-values (no correction for multiple testing)
#     pval <- sapply(1:ncol(soma_1), function(x) t.test(as.vector(unlist(soma_1[which(meta_1$Group == grp_1),x])), as.vector(unlist(soma_2[which(meta_2$Group == grp_2),x])) )$p.value)
#     
#     # correct for multiple testing
#     pval_adj <- p.adjust(pval, method = mtc)
#     alpha_adj <- alpha
#     if (mtc == "bonferroni"){
#       alpha_adj <- alpha / dim(soma)[2]
#     }
#     # still missing a method to plot some form of an adjusted significance level
#     # for BH correction
#     
#     # generate data frame
#     volc <- data.frame(log2fc = fc,
#                        pval = pval,
#                        pval_adj = pval_adj)
#     volc$difexpressed <- "NO"
#     volc$difexpressed[volc$log2fc > 0.6 & volc$pval_adj < alpha] <- "UP"
#     volc$difexpressed[volc$log2fc < -0.6 & volc$pval_adj < alpha] <- "DOWN"
#     volc$seqid = colnames(soma); rownames(volc) <- colnames(soma)
#     volc$gene <- sapply(1:ncol(soma), function(x) SomaScanAnnotation$Entrez.Gene.Name[SomaScanAnnotation$SeqId == colnames(soma)[x]])
#     n_de <- length(volc$difexpressed[!(volc$difexpressed == "NO")])
#     
#     plt_cols <- c("#00AFBB", "grey", "#bb0c00")
#     plt_labs <- c("Downregulated", "Not significant", "Upregulated")
#     if (length(unique(volc$difexpressed)) < 3){
#       if (sum(volc$difexpressed == "UP") == 0){
#         plt_cols <- plt_cols[-3]
#         plt_labs <- plt_labs[-3]
#       }
#       if (sum(volc$difexpressed == "DOWN") == 0){
#         plt_cols <- plt_cols[-1]
#         plt_labs <- plt_labs[-1]
#       }
#     }
#     
#     ## Plotting
#     # generate volcano plot
#     # theme_set(theme_classic(base_size = 12) +
#     #             theme(
#     #               axis.title.y = element_text(face = "bold", margin = margin(0,0,0,0), size = rel(1.1), color = 'black'),
#     #               axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(0,0,0,0), size = rel(1.1), color = 'black'),
#     #               plot.title = element_text(hjust = 0.5),
#     #               plot.margin = unit(c(0,0,0,0), 'lines')
#     #             ))
#     theme_set(theme_classic() +
#                 theme(plot.title = element_text(hjust = 0.5),
#                       legend.direction = legend.dir))
#     p <- volc %>%
#       ggplot(aes(x = log2fc,
#                  y = -log10(pval),
#                  col = difexpressed)) +
#       geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
#       geom_hline(yintercept = -log10(alpha_adj), col = "gray", linetype = 'dashed') +
#       geom_point(size = 2) +
#       scale_color_manual(values = plt_cols,
#                          labels = plt_labs) +
#       labs(color = 'Significance', #legend_title,
#            x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
#       geom_text_repel(data = subset(volc, !(difexpressed == "NO")),
#                       aes(label = gene)) +
#       scale_x_continuous(breaks = seq(-round(max(abs(c(min(volc$log2fc), max(volc$log2fc)))),1), round(max(abs(c(min(volc$log2fc), max(volc$log2fc)))),1), round(round(max(abs(c(min(volc$log2fc), max(volc$log2fc)))),1)/4,1))) +
#       annotation_custom(grobTree(textGrob(paste0("Correction: ", mtc_title, ",\n",
#                                                  "DE targets: ", n_de, ",\n",
#                                                  "Alpha = ", alpha), x = 0.02, y = 0.93, hjust = 0,
#                                           gp = gpar(col = "black")))) +
#     ggtitle(paste0(unique(meta_1$SampleGroup[meta_1$Group == grp_1])))
#     
#     # Store plot
#     myplots[[i]] <- local(p)
#   }
#   
#   legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
#   arranged_plots <- ggarrange(plotlist = myplots,
#                               ncol = length(grp_1_list),
#                               common.legend = TRUE,
#                               legend.grob = legend,
#                               hjust = -2)
#   title <- paste0(datanames[1], " vs ", datanames[2])
#   if (norm == TRUE){
#     title <- paste0(title, ", cross-normalised")
#   } else if (norm == "limma"){
#     title <- paste0(title, ", batch effects corrected (limma)")
#   }
#   annotate_figure(arranged_plots,
#                   top = text_grob(title,
#                                   face = "bold", size = 14))
#   
#   filename <- paste0("figures/", datanames[1], "_", datanames[2], "_volcano")
#   if (norm == TRUE){
#     filename <- paste0(filename, "_crossnorm")
#   }
#   if (norm == "limma"){
#     filename <- paste0(filename, "_limma")
#   }
#   if (mtc != "none"){
#     filename <- paste0(filename, "_", mtc)
#   }
#   if (alpha != 0.05){
#     filename <- paste0(filename, "_alph_", alpha)
#   }
#   filename <- paste0(filename, ".png")
#   ggsave(filename, width = 4000, height = 900, units = "px", dpi = 200, bg = "white")
# }


# ################################################################################
# ## Looking at calibrator-based cross-normalisation
# legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
# for (lesgo in "go"){
#   # Plot
#   myplots <- list(volcano_plot(soma1 = cal_norm_1_4_JYNR210602, soma2 = log_1_4_JYNR220601, meta2 = meta_mouse_JYNR220601,
#                                group1 = 1, group2 = 1, title = "Chow Vehicle"),
#                   volcano_plot(soma1 = cal_norm_1_4_JYNR210602, soma2 = log_1_4_JYNR220601, meta2 = meta_mouse_JYNR220601,
#                                group1 = 3, group2 = 2, title = "DIO-NASH Vehicle"))
#   
#   g <- ggarrange(plotlist = myplots,
#                  ncol = length(myplots),
#                  common.legend = TRUE,
#                  legend.grob = legend,
#                  hjust = -2)
#   title <- "JYNR210602 vs. JYNR220601, calibrator cross-normalised"
#   title <- paste0(title, ", one scale factor")
#   # title <- paste0(title, ", three scale factors")
#   # title <- paste0(title, ", individual scale factors")
#   title <- paste0(title, ", 1-4 steps only")
#   annotate_figure(g, top = text_grob(title))
#   filename <- "figures/JYNR210602_JYNR220601_volcano_contr_calnorm_"
#   filename <- paste0(filename, "one", "_1_4_bh.png")
#   # filename <- paste0(filename, "three", "_1_4_bh.png")
#   # filename <- paste0(filename, "indi", "_1_4_bh.png")
#   ggsave(filename, width = 4000, height = 1200, units = "px", dpi = 200, bg = "white")
# }



# ################################################################################
# ## Checking 1-4 vs. 1-5 and fc_threshold differences
# 
# for (lesgo in "go"){
#   legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
#   plotlist <- list()
#   for (fc in c(0, 0.3, 0.6, 1)){
#     for (handle in c("log_", "log_1_5_")){
#       soma_1 <- get(paste0(handle, "JYNR210602"))
#       soma_2 <- get(paste0(handle, "JYNR220601"))
#       if (handle == "log_"){
#         norm <- "partially normalised"
#       } else {
#         norm <- "fully normalised"
#       }
#       title_fc <- paste0(", log2(FC) threshold = ", fc)
#       myplots <- list(volcano_plot(soma1 = soma_1, title = paste0("JYNR210602, ", norm, title_fc), fc_cutoff = fc),
#                       volcano_plot(soma1 = soma_2, meta1 = meta_mouse_JYNR220601, group1 = 2, title = paste0("JYNR220601, ", norm, title_fc), fc_cutoff = fc))
#       plotlist <- append(plotlist, myplots)
#     }
#     
#     ggarrange(plotlist = plotlist,
#               ncol = length(myplots),
#               nrow = length(myplots),
#               common.legend = TRUE,
#               legend.grob = legend,
#               hjust = -2)
#     # title <- "JYNR210602 (top) vs. JYNR220602 (bottom), steps 1-4 (left) vs. steps 1-5 (right)"
#     # annotate_figure(g, top = text_grob(title))
#     filename <- "figures/JYNR210602_JYNR220601_volc_1_4v1_5_fc_"
#     filename <- paste0(filename, fc)
#     filename <- paste0(filename, ".png")
#     ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 185, bg = "white")
#   }
# }



################################################################################
## Z-score volcano plot
# SCRAP. MAKES NO SENSE TO DO.
# for (lesgo in "go"){
#   
#   log <- FALSE
#   log <- TRUE
#   
#   legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
#   soma_1 <- z_JYNR210602
#   soma_2 <- z_JYNR220601
#   
#   if (log == TRUE){
#     soma_1 <- z_log_JYNR210602; soma_2 <- z_log_JYNR220601
#   }
#   myplots <- list(volcano_plot(soma1 = soma_1, title = "JYNR210602, case vs. control"),
#                   volcano_plot(soma1 = soma_2, meta1 = meta_mouse_JYNR220601,
#                                title = "JYNR220601, case vs. control", group1 = 2),
#                   volcano_plot(soma1 = soma_1, soma2 = soma_2, meta2 = meta_mouse_JYNR220601,
#                                group1 = 1, group2 = 1, title = "Both studies, chow Vehicle"),
#                   volcano_plot(soma1 = soma_1, soma2 = soma_2, meta2 = meta_mouse_JYNR220601,
#                                group1 = 3, group2 = 2, title = "Both studies, DIO-NASH Vehicle"))
#   g <- ggarrange(plotlist = myplots,
#             ncol = 2, nrow = 2,
#             common.legend = TRUE,
#             legend.grob = legend,
#             hjust = -2)
#   title <- "Volcano plots based on Z-scores"
#   filename <- "figures/JYNR210602_JYNR220601_volc_z"
#   if (log == TRUE){
#     title <- paste0(title, " (log2-transformed)")
#     filename <- paste0(filename, "_log")
#   }
#   annotate_figure(g, top = text_grob(title))
#   filename <- paste0(filename, ".png")
#   ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 185, bg = "white")
# }





# ###############################################################################
# ## Outlier removal plot
# # Individual + inter-study plots
# for (lesgo in "go"){
# 
#   legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
#   soma_1 <- outrm_JYNR210602
#   soma_2 <- outrm_JYNR220601
# 
#   myplots <- list(volcano_plot(soma1 = soma_1, meta1 = outrm_meta_JYNR210602,
#                                title = "JYNR210602, case vs. control"),
#                   volcano_plot(soma1 = soma_2, meta1 = outrm_meta_JYNR220601,
#                                title = "JYNR220601, case vs. control", group1 = 2),
#                   volcano_plot(soma1 = soma_1, soma2 = soma_2,
#                                meta1 = outrm_meta_JYNR210602, meta2 = outrm_meta_JYNR220601,
#                                group1 = 1, group2 = 1, title = "Both studies, chow Vehicle"),
#                   volcano_plot(soma1 = soma_1, soma2 = soma_2, 
#                                meta1 = outrm_meta_JYNR210602, meta2 = outrm_meta_JYNR220601,
#                                group1 = 3, group2 = 2, title = "Both studies, DIO-NASH Vehicle"))
#   g <- ggarrange(plotlist = myplots,
#             ncol = 2, nrow = 2,
#             common.legend = TRUE,
#             legend.grob = legend,
#             hjust = -2)
#   title <- "Volcano plots with outliers removed"
#   filename <- "figures/JYNR210602_JYNR220601_volc_outrm"
#   annotate_figure(g, top = text_grob(title))
#   filename <- paste0(filename, ".png")
#   ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 185, bg = "white")
# }
# 
# # Plotting sample removal from JYNR210602, one by one
# for (go in "lesgo"){
#   grp1 <- meta_mouse_JYNR210602$Group == 3; grp2 <- meta_mouse_JYNR220601$Group == 1
#   myplots <- list(volcano_plot(soma1 = log_JYNR210602[grp1,][-1,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #1"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-2,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #2"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-3,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #3"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-4,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #4"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-5,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #5"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-6,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #6"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-7,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #7"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-8,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #8"))
#   # So apparently this one outlier pulled the entire grp. 3 of JYNR210602
#   # towards grp. 2 of JYNR220601 because removing it creates an abyss between them.
#   # Nah, it's just that its removal decreases the CV% and thus increases all p-values (higher power).
#   ggarrange(plotlist = myplots,
#             ncol = 4, nrow = 2,
#             common.legend = TRUE,
#             legend.grob = legend,
#             hjust = -2)
#   ggsave("figures/JYNR210602_JYNR220601_volc_outrm_1_8.png", width = 4000, height = 2000, units = "px", dpi = 185, bg = "white")
# }
# 
# # Individual + swap plots
# for (lesgo in "go"){
#   
#   legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
#   soma_1 <- outrm_JYNR210602
#   soma_2 <- outrm_JYNR220601
#   
#   myplots <- list(volcano_plot(soma1 = soma_1, meta1 = outrm_meta_JYNR210602,
#                                title = "JYNR210602, case vs. control"),
#                   volcano_plot(soma1 = soma_2, meta1 = outrm_meta_JYNR220601,
#                                title = "JYNR220601, case vs. control", group1 = 2),
#                   volcano_plot(soma1 = soma_1, soma2 = soma_2,
#                                meta1 = outrm_meta_JYNR210602, meta2 = outrm_meta_JYNR220601,
#                                group1 = 3, group2 = 1, title = "JYNR210602, placebo swapped"),
#                   volcano_plot(soma1 = soma_2, soma2 = soma_1, 
#                                meta1 = outrm_meta_JYNR220601, meta2 = outrm_meta_JYNR210602,
#                                group1 = 2, group2 = 1, title = "JYNR220601, placebo swapped"))
#   g <- ggarrange(plotlist = myplots,
#                  ncol = 2, nrow = 2,
#                  common.legend = TRUE,
#                  legend.grob = legend,
#                  hjust = -2)
#   title <- "Volcano plots with outliers removed"
#   filename <- "figures/JYNR210602_JYNR220601_volc_outrm_swap"
#   annotate_figure(g, top = text_grob(title))
#   filename <- paste0(filename, ".png")
#   ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 185, bg = "white")
# }







# ###############################################################################
# ## Outlier removal plot, version 2 (removed pre-normalisation)
# # Individual + inter-study plots
# for (lesgo in "go"){
#   
#   legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
#   soma_1 <- outrm2_JYNR210602
#   soma_2 <- outrm2_JYNR220601
#   
#   myplots <- list(volcano_plot(soma1 = soma_1, meta1 = outrm_meta_JYNR210602,
#                                title = "JYNR210602, case vs. control"),
#                   volcano_plot(soma1 = soma_2, meta1 = outrm_meta_JYNR220601,
#                                title = "JYNR220601, case vs. control", group1 = 2),
#                   volcano_plot(soma1 = soma_1, soma2 = soma_2,
#                                meta1 = outrm_meta_JYNR210602, meta2 = outrm_meta_JYNR220601,
#                                group1 = 1, group2 = 1, title = "Both studies, chow Vehicle"),
#                   volcano_plot(soma1 = soma_1, soma2 = soma_2, 
#                                meta1 = outrm_meta_JYNR210602, meta2 = outrm_meta_JYNR220601,
#                                group1 = 3, group2 = 2, title = "Both studies, DIO-NASH Vehicle"))
#   g <- ggarrange(plotlist = myplots,
#                  ncol = 2, nrow = 2,
#                  common.legend = TRUE,
#                  legend.grob = legend,
#                  hjust = -2)
#   title <- "Volcano plots with outliers removed (pre-normalisation)"
#   filename <- "figures/JYNR210602_JYNR220601_volc_outrm2"
#   annotate_figure(g, top = text_grob(title))
#   filename <- paste0(filename, ".png")
#   ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 185, bg = "white")
# }
# 
# # Plotting sample removal from JYNR210602, one by one
# for (go in "lesgo"){
#   myplots <- list(volcano_plot(soma1 = log_JYNR210602[grp1,][-1,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #1"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-2,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #2"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-3,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #3"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-4,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #4"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-5,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #5"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-6,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #6"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-7,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #7"),
#                   volcano_plot(soma1 = log_JYNR210602[grp1,][-8,], soma2 = log_JYNR220601[grp2,],
#                                meta_aligned = TRUE, title = "Missing sample #8"))
#   # So apparently this one outlier pulled the entire grp. 3 of JYNR210602
#   # towards grp. 2 of JYNR220601 because removing it creates an abyss between them.
#   ggarrange(plotlist = myplots,
#             ncol = 4, nrow = 2,
#             common.legend = TRUE,
#             legend.grob = legend,
#             hjust = -2)
#   ggsave("figures/JYNR210602_JYNR220601_volc_outrm2_1_8.png", width = 4000, height = 2000, units = "px", dpi = 185, bg = "white")
# }
# 
# # Individual + swap plots
# for (lesgo in "go"){
#   
#   legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
#   soma_1 <- outrm2_JYNR210602
#   soma_2 <- outrm2_JYNR220601
#   
#   myplots <- list(volcano_plot(soma1 = soma_1, meta1 = outrm_meta_JYNR210602,
#                                title = "JYNR210602, case vs. control"),
#                   volcano_plot(soma1 = soma_2, meta1 = outrm_meta_JYNR220601,
#                                title = "JYNR220601, case vs. control", group1 = 2),
#                   volcano_plot(soma1 = soma_1, soma2 = soma_2,
#                                meta1 = outrm_meta_JYNR210602, meta2 = outrm_meta_JYNR220601,
#                                group1 = 3, group2 = 1, title = "JYNR210602, placebo swapped"),
#                   volcano_plot(soma1 = soma_2, soma2 = soma_1, 
#                                meta1 = outrm_meta_JYNR220601, meta2 = outrm_meta_JYNR210602,
#                                group1 = 2, group2 = 1, title = "JYNR220601, placebo swapped"))
#   g <- ggarrange(plotlist = myplots,
#                  ncol = 2, nrow = 2,
#                  common.legend = TRUE,
#                  legend.grob = legend,
#                  hjust = -2)
#   title <- "Volcano plots with outliers removed (pre-normalisation)"
#   filename <- "figures/JYNR210602_JYNR220601_volc_outrm2_swap"
#   annotate_figure(g, top = text_grob(title))
#   filename <- paste0(filename, ".png")
#   ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 185, bg = "white")
# }
