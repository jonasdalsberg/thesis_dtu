## Initialize document
setwd("~/JDJG/clean")
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel)
source("scripts/SomaScan_Mouse/overlap_function.R")
source("scripts/SomaScan_Mouse/annotation_compass_function.R")

## Load data and perform initial wrangling (this should be further automated for future use)
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")

# Designed to plot confusion matrices, based on a reference list
for (letsgo in "plot"){
  datasets <- c("JYNR210602", "JYNR220601")
  ref_somas <- list(log_JYNR210602, log_JYNR220601)
  metas <- list(meta_mouse_JYNR210602, meta_mouse_JYNR220601)
  cases <- c(3,2)

  one_to_four <- TRUE
  alpha <- 0.05
  fc_cutoff <- 0
  fc_cutoff <- 0.3
  
  for (k in 1:2){
    # for (norm in c(FALSE, TRUE, "limma")){
    for (norm in c(FALSE)){
      # Parameters
      dataset1 <- datasets[k]; dataset2 <- rev(datasets)[k]
      ref_soma <- ref_somas[[k]];
      ref_meta <- meta1 <- metas[[k]]; meta2 <- rev(metas)[[k]];
      case1 <- cases[k]; contr <- 1


      title <- paste0(dataset1, ", confusion matrix")
      if (norm == TRUE){
        soma1 <- "crossnorm_JYNR210602"; soma2 <- "crossnorm_JYNR220601"
        title <- paste0(title, ", cross-normalised")
        print("cross-norm")
      } else if (norm == "limma"){
        soma1 <- "limma_JYNR210602"; soma2 <- "limma_JYNR220601"
        title <- paste0(title, ", batch effects corrected (limma)")
        print ("limma")
        if (one_to_four == TRUE){
          soma1 <- "limma_JYNR210602"; soma2 <- "limma_JYNR220601"
          print("1-4 steps only")
        }
      } else {
        soma1 <- "log_JYNR210602"; soma2 <- "log_JYNR220601"
        print("log")
        if (one_to_four == TRUE){
          soma1 <- "log_JYNR210602"; soma2 <- "log_JYNR220601"
          print("1-4 steps only")
        }
      }
      soma1 <- get(paste0(scaled, soma1)); soma2 <- get(paste0(scaled, soma2))

      if (one_to_four == TRUE){
        title <- paste0(title, ", 1-4 steps only")
      }
      if (k == 2){
        somas <- list(soma1, soma2)
        soma1 <- somas[[2]]; soma2 <- somas[[1]]
      }


      p <- plot_overlaps(soma1 = soma1, meta1 = meta1, case1 = case1, control1 = contr,
                         soma2 = soma2, meta2 = meta2, control2 = contr,
                         ref_soma = ref_soma, ref_meta = ref_meta,
                         ref_case = case1, ref_contr = contr,
                         alpha = alpha, mtc = mtc, title = title, fc_cutoff = fc_cutoff)

      p
      filename <- paste0("figures/", dataset1, "_", dataset2, "_cm")
      if (norm == TRUE){
        filename <- paste0(filename, "_crossnorm")
      }
      if (norm == "limma"){
        filename <- paste0(filename, "_limma")
      }
      # if (scaled == "scaled_"){
      #   filename <- paste0(filename, "_scaled")
      # }
      if (one_to_four == TRUE){
        filename <- paste0(filename, "")
      }
      if (!mtc == "none"){
        filename <- paste0(filename, "_", mtc)
      }
      if (alpha != 0.05){
        filename <- paste0(filename, "_alph_", alpha)
      }
      if (fc_cutoff != 0.3){
        filename <- paste0(filename, "_fc_", fc_cutoff)
      }
      filename <- paste0(filename, ".png")
      ggsave(filename, width = 2000, height = 2400, units = "px", dpi = 200, bg = "white")
    }
  }
}

for (df in test){
  test2 <- df[1]
}

# Designed to plot confusion matrices
# for (norm in c(FALSE, TRUE, "limma")){
for (norm in c(FALSE)){
  for (datasets in list(c("JYNR210602", "JYNR220601"),
                        c("JYNR220601", "JYNR210602"))){
    for (one_to_four in c(TRUE, FALSE)){
      # Parameters
      dataset1 <- datasets[1]; dataset2 <- datasets[2]
      # norm <- TRUE
      # norm <- "limma"
      # norm <- FALSE
      alpha <- 0.05
      # alpha <- 0.1
      mtc <- "none"
      mtc <- "bh"
      fc_cutoff <- 0
      fc_cutoff <- 0.3
      
      
      
      case1 <- 3; control1 <- 1; case2 <- 2; control2 <- 1
      meta1 <- meta_mouse_JYNR210602; meta2 <- meta_mouse_JYNR220601
      title <- paste0(dataset1, ", confusion matrix")
      if (norm == TRUE){
        soma1 <- "crossnorm_JYNR210602"; soma2 <- "crossnorm_JYNR220601"
        title <- paste0(title, ", cross-normalised")
        print("cross-norm")
      } else if (norm == "limma"){
        soma1 <- "limma_JYNR210602"; soma2 <- "limma_JYNR220601"
        title <- paste0(title, ", batch effects corrected (limma)")
        print ("limma")
        if (one_to_four == TRUE){
          soma1 <- "limma_JYNR210602"; soma2 <- "limma_JYNR220601"
          print("1-4 steps only")
        }
      } else {
        soma1 <- "log_1_5_JYNR210602"; soma2 <- "log_1_5_JYNR220601"
        print("log")
        if (one_to_four == TRUE){
          soma1 <- "log_JYNR210602"; soma2 <- "log_JYNR220601"
          print("1-4 steps only")
        }
      }
      soma1 <- get(paste0(scaled, soma1)); soma2 <- get(paste0(scaled, soma2))
      
      if (one_to_four == TRUE){
        title <- paste0(title, ", 1-4 steps only")
      }
      if (fc_cutoff != 0.3){
        title <- paste0(title, ", log2(FC) threshold = ", fc_cutoff)
      }
      
      if (dataset1 == "JYNR220601"){
        soma_temp <- soma1; soma1 <- soma2; soma2 <- soma_temp
        meta_temp <- meta1; meta1 <- meta2; meta2 <- meta_temp
        case_temp <- case1; case1 <- case2; case2 <- case_temp
      }
      
      p <- plot_overlaps(soma1 = soma1, meta1 = meta1, case1 = case1, control1 = control1,
                         soma2 = soma2, meta2 = meta2, control2 = control2,
                         alpha = alpha, mtc = mtc, title = title, fc_cutoff = fc_cutoff)
      
      p
      filename <- paste0("figures/", dataset1, "_", dataset2, "_cm")
      if (norm == TRUE){
        filename <- paste0(filename, "_crossnorm")
      }
      if (norm == "limma"){
        filename <- paste0(filename, "_limma")
      }
      # if (scaled == "scaled_"){
      #   filename <- paste0(filename, "_scaled")
      # }
      if (one_to_four == TRUE){
        filename <- paste0(filename, "_1_4")
      }
      if (!mtc == "none"){
        filename <- paste0(filename, "_", mtc)
      }
      if (alpha != 0.05){
        filename <- paste0(filename, "_alph_", alpha)
      }
      if (fc_cutoff != 0.3){
        filename <- paste0(filename, "_fc_", fc_cutoff)
      }
      filename <- paste0(filename, ".png")
      ggsave(filename, width = 2000, height = 2400, units = "px", dpi = 200, bg = "white")
      
    }
  }
}


# ################################################################################
# # This chunk is designed to plot confusion matrices for the realistic normalisation dfs
# for (get_all_case_control_norm_placebo_swapped_plots in "go"){
#   # Plot and save (designed for plotting the best contrasts in both studies)
#   legend <- get_legend(volcano_plot(legend.dir = "horizontal"))
#   n_plots <- length(case_control_norm_dfs)/2 / 2 # skips scenarios 3-4
#   norms <- c(TRUE, "limma", TRUE, "limma", FALSE, "limma", FALSE, "limma")
#   steps_incl <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
#   scenario <- c(1,1,2,2,3,3,4,4)
#   data_cases <- c("JYNR210602", "JYNR220601")
#   data_controls <- c("JYNR220601", "JYNR210602")
#   ref_somas <- list(log_JYNR210602, log_JYNR220601)
#   ref_metas <- list(meta_mouse_JYNR210602, meta_mouse_JYNR220601)
#   ref_cases <- c(3,2)
#   # dataset2 <- "JYNR210602"; dataset1 <- "JYNR220601" # comment out to plot JYNR210602
#   for (j in 1:n_plots){
#     for (k in 1:2){
#       print(j)
#       i <- j*2 - 2 + k # gives numbers 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16
#       # Parameters
#       # norm <- TRUE
#       # norm <- "limma"
#       # norm <- FALSE
#       alpha <- 0.05
#       # alpha <- 0.1
#       mtc <- "none"
#       mtc <- "bh"
#       
#       # Choose data
#       case <- case_control_norm_dfs[[i]][[1]]; placebo <- case_control_norm_dfs[[i]][[2]]
#       ref_soma <- ref_somas[[k]]; ref_meta <- ref_metas[[k]]; ref_case <- ref_cases[k]
#       norm <- norms[j]; one_to_four <- steps_incl[j]
#       print("Data chosen.")
#       
#       # Dataset 1 (JYRN210602)
#       # Plot confusion matrices
#       dataset1 <- data_cases[k]; dataset2 <- data_controls[k]
#       title <- paste0(dataset1, ", confusion matrix")
#       if (norm == TRUE){
#         title <- paste0(title, ", cross-normalised")
#       }
#       if (norm == "limma"){
#         title <- paste0(title, ", batch effects corrected (limma)")
#       }
#       if (one_to_four == TRUE){
#         title <- paste0(title, ", 1-4 steps only")
#       }
#       title <- paste0(title, ", scenario ", scenario[j])
#       p <- plot_overlaps(soma1 = case, soma2 = placebo, meta_aligned = TRUE,
#                          ref_soma = ref_soma, ref_meta = ref_meta,
#                          ref_case = ref_case, ref_contr = 1,
#                          alpha = alpha, mtc = mtc, title = title)
#       
#       p
#       filename <- paste0("figures/", dataset1, "_", dataset2, "_cm")
#       if (norm == TRUE){
#         filename <- paste0(filename, "_crossnorm")
#       }
#       if (norm == "limma"){
#         filename <- paste0(filename, "_limma")
#       }
#       if (one_to_four == TRUE){
#         filename <- paste0(filename, "_1_4")
#       }
#       if (!mtc == "none"){
#         filename <- paste0(filename, "_", mtc)
#       }
#       if (alpha != 0.05){
#         filename <- paste0(filename, "_alph_", alpha)
#       }
#       filename <- (paste0(filename, "_sce", scenario[j]))
#       filename <- paste0(filename, ".png")
#       ggsave(filename, width = 2000, height = 2400, units = "px", dpi = 200, bg = "white")
#       print(paste0("Confusion matrix saved (dataset 1), scenario ", scenario[j], ", norm = ", norm, ", 1-4 = ", one_to_four)) 
#     }
#   }
# }




################################################################################
# Determine optimal conditions (designed to optimise precision)
# Uses JYNR210602 (case/control) with JYNR220601 (control swap)
for (run_this_line_to_run_everything in "yes"){
  alphas <- c(0.05)
  # alphas <- c(0.05, 10^-2, 10^-3, 10^-4)
  # alphas <- c(0.5, 0.2, 0.1, 0.05, 10^-5)
  fc_cutoffs <- c(0.3)
  fc_cutoffs <- c(0.3, 0.6, 0.9, 1.2)
  # fc_cutoffs <- c(0, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 4)
  # norms <- c("log_", "crossnorm_", "limma_")
  norms <- c("log_")
  combinations <- length(fc_cutoffs) * length(alphas) * length(norms) * 2
  store_data <- data.frame(dataset = c(rep("JYNR210602", combinations/2), rep("JYNR220601", combinations/2)),
                           norm = rep(NA, combinations),
                           FC_cutoff = rep(NA, combinations),
                           alpha = rep(NA, combinations),
                           precision = rep(NA, combinations),
                           positives = rep(NA, combinations))
  cases <- c("JYNR210602", "JYNR220601"); controls <- rev(cases)
  metas <- list(meta_mouse_JYNR210602, meta_mouse_JYNR220601)
  group1 <- c(3,2); contr <- 1

  refs <- list(get_signif_targets(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                                  group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = 0.3),
               get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                                  group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = 0.3))

  i <- 0
  for (j in 1:2){
    dataset <- cases[j]; print(dataset)
    ref <- refs[[j]]
    case <- group1[j]
    meta1 <- metas[[j]]; meta2 <- rev(metas)[[j]]
    for (norm in norms){
      # Define datasets
      if (norm == "norm_"){
        normalisation <- "Cross-norm"
        print("Cross-normalization (step 5)")
      } else if (norm == "limma_"){
        normalisation <- "Limma"
        print ("Batch effect correction (Limma)")
      } else {
        normalisation <- "Baseline"
        print("Baseline")
      }
      soma1 <- get(paste0(norm, cases[j])); soma2 <- get(paste0(norm, controls[j]))

      for(alpha in alphas){
        for (fc_cutoff in fc_cutoffs){
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

          store_data[i,] <- c(dataset,
                              normalisation,
                              fc_cutoff,
                              alpha,
                              precision,
                              length(positive))
        }
      }
    }
  }
}

## Plot parameters
# Precision vs. alpha
# alpha_fc0_data <- store_data
# alpha_store_data <- store_data
alpha_data <- alpha_store_data
alpha_data[,4:6] <- lapply(4:6, function(x) as.numeric(alpha_data[,x]))
# ymin <- min(alpha_data$precision)*0.85; ymax <- max(alpha_data$precision)
a1 <- alpha_data %>% ggplot(mapping = aes(x = alpha,
                                    y = precision,
                                    colour = dataset)) +
  ylim(c(0,1)) +
  geom_line() +
  ggtitle("Significant proteins after placebo swapping")+
  scale_x_reverse() +
  xlab("Significance threshold (alpha)") +
  ylab("Precision (TP/(TP+FP))") +
  labs(colour = "Dataset") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))   +
  annotation_compass("log2(FC) threshold = 0.3", "SW")

# TP vs. alpha
ymax1 <- max(alpha_data$positives)
a2 <- alpha_data %>% ggplot(mapping = aes(x = alpha,
                                          y = positives,
                                          color = dataset)) +
  geom_line() +
  ggtitle("Significant proteins after placebo swapping") +
  scale_x_reverse() +
  xlab("Significance threshold (alpha)") +
  ylab("DE targets after swapping") +
  labs(colour = "Dataset") +
  theme_classic() +
  ylim(0,ymax1) +
  theme(plot.title = element_text(hjust = 0.5))  +
  annotation_compass("log2(FC) threshold = 0.3", "SW")

# ggarrange(a1, a2)

# Precision vs. FC_cutoff
# fc_store_data <- store_data
fc_data <- fc_store_data
# fc_data$precision[45] <- 0
fc_data[,3:6] <- lapply(3:6, function(x) as.numeric(fc_data[,x]))
# ymin <- min(fc_data$precision)*0.85; ymax <- max(fc_data$precision)
f1 <- fc_data %>% ggplot(mapping = aes(x = FC_cutoff,
                                       y = precision,
                                       color = dataset)) +
  ylim(c(0,1)) +
  geom_line() +
  ggtitle("Significant proteins after placebo swapping") +
  xlab("Log2(fold change) threshold") +
  ylab("Precision (TP/(TP+FP))") +
  labs(colour = "Dataset") +
  ylim(c(0,1)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotation_compass("alpha = 0.05", "SW")

# TP vs. alpha
ymax2 <- max(fc_data$positives)
x1_start <- min(fc_data$FC_cutoff); x2_end <- max(fc_data$FC_cutoff)
f2 <- fc_data %>% ggplot(mapping = aes(x = FC_cutoff,
                                       y = positives,
                                       color = dataset)) +
  geom_line() +
  ggtitle("Significant proteins after placebo swapping") +
  xlab("Log2(fold change) threshold") +
  ylab("Significant targets after swapping") +
  labs(colour = "Dataset") +
  theme_classic() +
  ylim(c(0,ymax2)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotation_compass("alpha = 0.05", "NE") +
  geom_segment(x = x1_start, xend = x2_end, y = 100, yend = 100,
               linetype = "dashed", colour = "lightgrey") +
  geom_segment(x = x1_start, xend = x2_end, y = 50, yend = 50,
               linetype = "dashed", colour = "darkgrey") +
  geom_segment(x = x1_start, xend = x2_end, y = 0, yend = 0,
               linetype = "dashed", colour = "black")

ggarrange(a1, f1, a2, f2)

ggsave("figures/optimise_alpha_fc_2.png", width = 800*3, height = 420*3, units = "px", dpi = 200, bg = "white")



################################################################################
# Designed to plot precision for different normalisation methods
for (plot in "letsgo"){
  # Define parameters
  fc_cutoff = 0
  fc_cutoff = 0.3
  test_dif_qualities = FALSE
  test_dif_qualities = TRUE
  
  
  
  if (test_dif_qualities == TRUE){
    quals <- 1:3
    seq_qual <- list(
      SomaScanAnnotation$SeqId[SomaScanAnnotation$mouse_quality %in% c("Low") & SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)],
      SomaScanAnnotation$SeqId[SomaScanAnnotation$mouse_quality %in% c("Medium_biosignal", "Medium_dilution") & SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)],
      SomaScanAnnotation$SeqId[SomaScanAnnotation$mouse_quality %in% c("High") & SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)]
    )
    qual_titles <- c("low quality aptamers", "medium quality aptamers", "high quality aptamers")
    qual_fileexs <- c("low", "med", "high")
  } else {
    quals <- 1
    seq_qual <- list(SomaScanAnnotation$SeqId[SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)])
  }
  
  for (q_i in quals){
    
    seqids <- seq_qual[[q_i]]
    
    cases <- list(list(
        log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # No normalisation (1-4)
        crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Cross-normalisation
        case_control_norm_dfs[[1]][[1]], # scenario 1, cross-normalised (JYNR210602)
        case_control_norm_dfs[[5]][[1]], # scenario 2, cross-normalised (JYNR210602)
        limma_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Limma (1-4)
        case_control_norm_dfs[[7]][[1]], # scenario 2, limma (JYNR210602)
        cal_one_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Calibrator-normalised (1 SF)
        cal_three_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Calibrator-normalised (3 SFs)
        cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 3,]), # Calibrator-normalised (indi SFs)
      
      list(
        log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
        crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Cross-normalisation
        case_control_norm_dfs[[2]][[1]], # scenario 1, cross-normalised (JYNR220601)
        case_control_norm_dfs[[6]][[1]], # scenario 2, cross-normalised (JYNR220601)
        limma_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Limma (1-4)
        case_control_norm_dfs[[8]][[1]], # scenario 2, limma (JYNR220601)
        log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Calibrator-normalised (1 SF)
        log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Calibrator-normalised (3 SFs)
        log_JYNR220601[meta_mouse_JYNR220601$Group == 2,] #  Calibrator-normalised (indi SFs)
    ))
    
    placebos <- list(list(
        log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
        crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Cross-normalisation
        case_control_norm_dfs[[1]][[2]], # scenario 1, cross-normalised (JYNR210602)
        case_control_norm_dfs[[5]][[2]], # scenario 5, cross-normalised (JYNR210602)
        limma_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Limma (1-4)
        case_control_norm_dfs[[7]][[2]], # scenario 5, limma (JYNR210602)
        log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Calibrator-normalised (1 SF)
        log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Calibrator-normalised (3 SFs)
        log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]), #  Calibrator-normalised (indi SFs)
      
      list(
        log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # No normalisation (1-4)
        crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Cross-normalisation
        case_control_norm_dfs[[2]][[2]], # scenario 1, cross-normalised (JYNR220601)
        case_control_norm_dfs[[6]][[2]], # scenario 5, cross-normalised (JYNR220601)
        limma_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Limma (1-4)
        case_control_norm_dfs[[8]][[2]], # scenario 5, limma (JYNR220601)
        cal_one_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Calibrator-normalised (1 SF)
        cal_three_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Calibrator-normalised (3 SFs)
        cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 1,] # Calibrator-normalised (indi SFs)
      ))
    
    # Compute references (intra-study significant proteins)
      # Four references: Two for each study, one with intra-study normalisation (1-5) and one without (1-4)
    refs <- list(get_signif_targets(soma1 = log_JYNR210602[,colnames(log_JYNR210602) %in% seqids], meta1 = meta_mouse_JYNR210602,
                                  group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
                 get_signif_targets(soma1 = log_JYNR220601[,colnames(log_JYNR220601) %in% seqids], meta1 = meta_mouse_JYNR220601,
                                    group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff))
    
    
    precisions <- TPs <- positives <- matrix(nrow = 2, ncol = length(cases[[1]]))
    precisions <- TPs <- positives <- as.data.frame(precisions)
    rownames(precisions) <- rownames(TPs) <- rownames(positives) <- datasets <- c("JYNR210602", "JYNR220601")
    colnames(precisions) <- colnames(TPs) <- colnames(positives) <- c("Baseline",
                                                                      "Cross",
                                                                      "Cross, real 1",
                                                                      "Cross, real 2",
                                                                      "Limma",
                                                                      "Limma, real",
                                                                      "Cal_one",
                                                                      "Cal_three",
                                                                      "Cal_all")
    
    for (k in 1:2){ # for each dataset (JYNR210602 and JYNR220601)
      for (i in 1:dim(precisions)[2]){
        case <- cases[[k]][[i]]; placebo <- placebos[[k]][[i]]
        # Filter for SOMAmer quality (if true)
        case <- case[,colnames(case) %in% seqids]; placebo <- placebo[,colnames(placebo) %in% seqids]
        ref <- refs[[k]]
        positive <- get_signif_targets(soma1 = case, soma2 = placebo, fc_cutoff = fc_cutoff,
                                        meta_aligned = TRUE, sig_list_only = TRUE)
      TP <- sum(ref %in% positive)
      FP <- sum(!positive %in% ref)
      precision <- TP/(TP+FP)
      precisions[k,i] <- precision
      TPs[k,i] <- TP
      positives[k,i] <- length(positive)
      }
    }
    
    # Compute optimal precision (overlap btw. refs, no swapping)
    ref1 <- refs[[1]]; ref2 <- refs[[2]]
    pre_ceiling1 <- sum(ref1 %in% ref2) / length(ref2)
    pre_ceiling2 <- sum(ref2 %in% ref1) / length(ref1)
    # precision may not be a good metric for computing ceilings for datasets of different sizes
    
    # Put into long format and join dfs
    precisions$Dataset <- TPs$Dataset <- positives$Dataset <- rownames(precisions)
    
    precisions_long <- precisions %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "Precision",
                   cols = colnames(precisions[,!colnames(precisions) == "Dataset"]))
    TPs_long <- TPs %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "TPs",
                   cols = colnames(positives[,!colnames(positives) == "Dataset"]))
    positives_long <- positives %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "Positives",
                   cols = colnames(positives[,!colnames(positives) == "Dataset"]))
    norm_data <- precisions_long %>% 
      inner_join(TPs_long) %>% 
      inner_join(positives_long)
    norm_data$`Normalisation method` <- factor(norm_data$`Normalisation method`,
                                                  levels = c("Baseline",
                                                             "Intra-study",
                                                             "Cross",
                                                             "Cross, real 1",
                                                             "Cross, real 2",
                                                             "Cal_one",
                                                             "Cal_three",
                                                             "Cal_all",
                                                             "Limma",
                                                             "Limma + intra-study",
                                                             "Limma, real"))
    
    # levels(norm_data$`Normalisation method`)
    
    
    p1 <- norm_data %>% 
      ggplot(mapping = aes(x = `Normalisation method`,
                           y = Precision,
                           fill = Dataset)) +
      geom_col(position = "dodge",
               width = 0.4) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1),
            legend.direction = "horizontal") +
      annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
               x = 0.6, y = 1, hjust = 0) +
      geom_hline(yintercept = c(pre_ceiling1, pre_ceiling2), linetype = "dashed", color = c("#F8766D", "#00BFC4"))
  
    p2 <- norm_data %>%
      ggplot(mapping = aes(x = `Normalisation method`,
                           y = Positives,
                           fill = Dataset)) +
      geom_col(position = "dodge",
               width = 0.4) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1)) +
      annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
               x = 9.4, y = max(c(length(ref1), length(ref2), max(norm_data$Positives))), hjust = 1) +
      ylab("Significant targets after swapping") +
      geom_hline(yintercept = c(length(ref1), length(ref2)), linetype = "dashed", color = c("#F8766D", "#00BFC4"))
    
    legend <- get_legend(p1)
    g <- ggarrange(p1,p2,
              ncol = 2,
              common.legend = TRUE,
              legend.grob = legend,
              hjust = -2)
    title <- "Evaluation of normalisation method feasibility"
    if (test_dif_qualities == TRUE){
      title <- paste0(title, ", ", qual_titles[q_i])
    }
    annotate_figure(g, top = text_grob(title))
    # Noter til mig selv:
    # Sæt længde/bredde, så der kan være to af disse på et PP.
    # Nr. 2 skal være med FC_cutoff = 0.3.
    # Vi skal lige have en mere realistisk version af Limma (se One Note)
      # ^ Tjek. Det har vi nu. Hvordan ser det ud nu?
    filename <- paste0("figures/norm_perf_fc_", fc_cutoff)
    if (test_dif_qualities == TRUE){
      filename <- paste0(filename, "_", qual_fileexs[q_i])
    }
    filename <- paste0(filename, ".png")
    ggsave(filename, width = 4000, height = 1200, units = "px", dpi = 200, bg = "white")
  }
}







################################################################################
# Designed to plot precision for different normalisation methods (version 2)
for (plot in "letsgo"){
  # Define parameters
  fc_cutoff <- 0
  fc_cutoff <- 0.3
  test_dif_qualities <- FALSE
  # test_dif_qualities <- TRUE
  small_setup <- FALSE
  # small_setup <- TRUE
  
  
  
  if (small_setup == TRUE){
    norm_methods <- c("Baseline",
                      "Cross",
                      "Limma")
  } else {
    norm_methods <- c("Baseline",
                      "Cross",
                      "Cross, real 1",
                      "Cross, real 2",
                      "Limma",
                      "Limma, real 1",
                      "Limma, real 2",
                      "Cal")
  }
  
  
  if (test_dif_qualities == TRUE){
    quals <- 1:3
    seq_qual <- list(
      SomaScanAnnotation$SeqId[SomaScanAnnotation$mouse_quality %in% c("Low") & SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)],
      SomaScanAnnotation$SeqId[SomaScanAnnotation$mouse_quality %in% c("Medium_biosignal", "Medium_dilution") & SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)],
      SomaScanAnnotation$SeqId[SomaScanAnnotation$mouse_quality %in% c("High") & SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)]
    )
    qual_titles <- c("low quality aptamers", "medium quality aptamers", "high quality aptamers")
    qual_fileexs <- c("low", "med", "high")
  } else {
    quals <- 1
    seq_qual <- list(SomaScanAnnotation$SeqId[SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)])
  }
  
  for (q_i in quals){
    
    seqids <- seq_qual[[q_i]]
    
    cases <- list(list(
      log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # No normalisation (1-4)
      crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Cross-normalisation
      case_control_norm_dfs[[1]][[1]], # scenario 1, cross-normalised (JYNR210602)
      case_control_norm_dfs[[5]][[1]], # scenario 2, cross-normalised (JYNR210602)
      limma_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Limma (1-4)
      case_control_norm_dfs[[3]][[1]], # scenario 1, limma (JYNR210602)
      case_control_norm_dfs[[7]][[1]], # scenario 2, limma (JYNR210602)
      cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 3,]), # Calibrator-normalised (indi SFs)
      
      list(
        log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
        crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Cross-normalisation
        case_control_norm_dfs[[2]][[1]], # scenario 1, cross-normalised (JYNR220601)
        case_control_norm_dfs[[6]][[1]], # scenario 2, cross-normalised (JYNR220601)
        limma_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Limma (1-4)
        case_control_norm_dfs[[4]][[1]], # scenario 1, limma (JYNR220601)
        case_control_norm_dfs[[8]][[1]], # scenario 2, limma (JYNR220601)
        log_JYNR220601[meta_mouse_JYNR220601$Group == 2,] #  Calibrator-normalised (indi SFs)
      ))
    
    placebos <- list(list(
      log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
      crossnorm_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Cross-normalisation
      case_control_norm_dfs[[1]][[2]], # scenario 1, cross-normalised (JYNR210602)
      case_control_norm_dfs[[5]][[2]], # scenario 2, cross-normalised (JYNR210602)
      limma_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Limma (1-4)
      case_control_norm_dfs[[3]][[2]], # scenario 1, limma (JYNR210602)
      case_control_norm_dfs[[7]][[2]], # scenario 2, limma (JYNR210602)
      log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]), #  Calibrator-normalised (indi SFs)
      
      list(
        log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # No normalisation (1-4)
        crossnorm_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Cross-normalisation
        case_control_norm_dfs[[2]][[2]], # scenario 1, cross-normalised (JYNR220601)
        case_control_norm_dfs[[6]][[2]], # scenario 2, cross-normalised (JYNR220601)
        limma_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Limma (1-4)
        case_control_norm_dfs[[4]][[2]], # scenario 2, limma (JYNR220601)
        case_control_norm_dfs[[8]][[2]], # scenario 2, limma (JYNR220601)
        cal_indi_JYNR210602[meta_mouse_JYNR210602$Group == 1,] # Calibrator-normalised (indi SFs)
      ))
    
    # Compute references (intra-study significant proteins)
    # Four references: Two for each study, one with intra-study normalisation (1-5) and one without (1-4)
    refs <- list(get_signif_targets(soma1 = log_JYNR210602[,colnames(log_JYNR210602) %in% seqids], meta1 = meta_mouse_JYNR210602,
                                    group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
                 get_signif_targets(soma1 = log_JYNR220601[,colnames(log_JYNR220601) %in% seqids], meta1 = meta_mouse_JYNR220601,
                                    group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff))
    
    if (small_setup == TRUE){
      cases[[1]] <- cases[[1]][c(1,2,5)]; cases[[2]] <- cases[[2]][c(1,2,5)]
      placebos[[1]] <- placebos[[1]][c(1,2,5)]; placebos[[2]] <- placebos[[2]][c(1,2,5)]
    }
    
    precisions <- TPs <- positives <- matrix(nrow = 2, ncol = length(cases[[1]]))
    precisions <- TPs <- positives <- as.data.frame(precisions)
    rownames(precisions) <- rownames(TPs) <- rownames(positives) <- datasets <- c("JYNR210602", "JYNR220601")
    colnames(precisions) <- colnames(TPs) <- colnames(positives) <- norm_methods
    
    for (k in 1:2){ # for each dataset (JYNR210602 and JYNR220601)
      for (i in 1:dim(precisions)[2]){
        case <- cases[[k]][[i]]; placebo <- placebos[[k]][[i]]
        # Filter for SOMAmer quality (if true)
        case <- case[,colnames(case) %in% seqids]; placebo <- placebo[,colnames(placebo) %in% seqids]
        ref <- refs[[k]]
        positive <- get_signif_targets(soma1 = case, soma2 = placebo, fc_cutoff = fc_cutoff,
                                       meta_aligned = TRUE, sig_list_only = TRUE)
        TP <- sum(ref %in% positive)
        FP <- sum(!positive %in% ref)
        precision <- TP/(TP+FP)
        precisions[k,i] <- precision
        TPs[k,i] <- TP
        positives[k,i] <- length(positive)
      }
    }
    
    # Compute optimal precision (overlap btw. refs, no swapping)
    ref1 <- refs[[1]]; ref2 <- refs[[2]]
    pre_ceiling1 <- sum(ref1 %in% ref2) / length(ref2)
    pre_ceiling2 <- sum(ref2 %in% ref1) / length(ref1)
    # precision may not be a good metric for computing ceilings for datasets of different sizes
    
    # Put into long format and join dfs
    precisions$Dataset <- TPs$Dataset <- positives$Dataset <- rownames(precisions)
    
    precisions_long <- precisions %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "Precision",
                   cols = colnames(precisions[,!colnames(precisions) == "Dataset"]))
    TPs_long <- TPs %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "TPs",
                   cols = colnames(positives[,!colnames(positives) == "Dataset"]))
    positives_long <- positives %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "Positives",
                   cols = colnames(positives[,!colnames(positives) == "Dataset"]))
    norm_data <- precisions_long %>% 
      inner_join(TPs_long) %>% 
      inner_join(positives_long)
    norm_data$`Normalisation method` <- factor(norm_data$`Normalisation method`,
                                               levels = norm_methods)
    
    # levels(norm_data$`Normalisation method`)
    
    
    p1 <- norm_data %>% 
      ggplot(mapping = aes(x = `Normalisation method`,
                           y = Precision,
                           fill = Dataset)) +
      geom_col(position = "dodge",
               width = 0.4) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1),
            legend.direction = "horizontal") +
      annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
               x = 0.6, y = 1, hjust = 0) #+
      # geom_hline(yintercept = c(pre_ceiling1, pre_ceiling2), linetype = "dashed", color = c("#F8766D", "#00BFC4"))
    
    p2 <- norm_data %>%
      ggplot(mapping = aes(x = `Normalisation method`,
                           y = Positives,
                           fill = Dataset)) +
      geom_col(position = "dodge",
               width = 0.4) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1)) +
      # annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
      #          x =  length(norm_methods) + 0.4, y = max(c(length(ref1), length(ref2), max(norm_data$Positives))), hjust = 1) +
      ylab("Significant targets after swapping") +
      geom_hline(yintercept = c(length(ref1), length(ref2)), linetype = "dashed", color = c("#F8766D", "#00BFC4"))
    
    legend <- get_legend(p1)
    g <- ggarrange(p1,p2,
                   ncol = 2,
                   common.legend = TRUE,
                   legend.grob = legend,
                   hjust = -2)
    title <- "Evaluation of normalisation method feasibility"
    if (test_dif_qualities == TRUE){
      title <- paste0(title, ", ", qual_titles[q_i])
    }
    annotate_figure(g, top = text_grob(title))
    # Noter til mig selv:
    # Sæt længde/bredde, så der kan være to af disse på et PP.
    # Nr. 2 skal være med FC_cutoff = 0.3.
    # Vi skal lige have en mere realistisk version af Limma (se One Note)
    # ^ Tjek. Det har vi nu. Hvordan ser det ud nu?
    filename <- paste0("figures/norm_perf_fc_", fc_cutoff)
    if (test_dif_qualities == TRUE){
      filename <- paste0(filename, "_", qual_fileexs[q_i])
    }
    if (small_setup == TRUE){
      filename <- paste0(filename, "_basic")
    } else {
      filename <- paste0(filename, "_complete")
    }
    filename <- paste0(filename, ".png")
    ggsave(filename, width = 4000, height = 1200, units = "px", dpi = 200, bg = "white")
  }
}












################################################################################
# Designed to plot precision for different normalisation methods (steps 1-4 vs 1-5)
# Individual truths
for (plot in "letsgo"){
  # Define parameters
  fc_cutoff = 0.3
    
  cases <- list(
    list(
      log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # No normalisation (1-4)
      log_1_5_JYNR210602[meta_mouse_JYNR210602$Group == 3,] # Intra-study normalisation (1-5)
      ),
    list(
      log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
      log_1_5_JYNR220601[meta_mouse_JYNR220601$Group == 2,] # Intra-study normalisation (1-5)
      ))
    
  placebos <- list(
      list(
        log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
        log_1_5_JYNR220601[meta_mouse_JYNR220601$Group == 1,] # Intra-study normalisation (1-5)
      ),
      list(
        log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # No normalisation (1-4)
        log_1_5_JYNR210602[meta_mouse_JYNR210602$Group == 1,] # Intra-study normalisation (1-5)
        ))
    
    # Compute references (intra-study significant proteins)
    # Four references: Two for each study, one with intra-study normalisation (1-5) and one without (1-4)
    refs <- list(
      list(
        get_signif_targets(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                           group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
        get_signif_targets(soma1 = log_1_5_JYNR210602, meta1 = meta_mouse_JYNR210602,
                           group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff)),
      list(
        get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                           group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
        get_signif_targets(soma1 = log_1_5_JYNR220601, meta1 = meta_mouse_JYNR220601,
                           group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff)))
    
    
    precisions <- TPs <- positives <- matrix(nrow = 2, ncol = length(cases[[1]]))
    precisions <- TPs <- positives <- as.data.frame(precisions)
    rownames(precisions) <- rownames(TPs) <- rownames(positives) <- datasets <- c("JYNR210602", "JYNR220601")
    colnames(precisions) <- colnames(TPs) <- colnames(positives) <- c("Baseline",
                                                                      "Intra-study")
    
    for (k in 1:2){ # for each dataset (JYNR210602 and JYNR220601)
      for (i in 1:dim(precisions)[2]){
        case <- cases[[k]][[i]]; placebo <- placebos[[k]][[i]]
        ref <- refs[[k]][[i]]
        positive <- get_signif_targets(soma1 = case, soma2 = placebo, fc_cutoff = fc_cutoff,
                                       meta_aligned = TRUE, sig_list_only = TRUE)
        TP <- sum(ref %in% positive)
        FP <- sum(!positive %in% ref)
        precision <- TP/(TP+FP)
        precisions[k,i] <- precision
        TPs[k,i] <- TP
        positives[k,i] <- length(positive)
      }
    }
    
    
    # Put into long format and join dfs
    precisions$Dataset <- TPs$Dataset <- positives$Dataset <- rownames(precisions)
    
    precisions_long <- precisions %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "Precision",
                   cols = colnames(precisions[,!colnames(precisions) == "Dataset"]))
    TPs_long <- TPs %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "TPs",
                   cols = colnames(positives[,!colnames(positives) == "Dataset"]))
    positives_long <- positives %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "Positives",
                   cols = colnames(positives[,!colnames(positives) == "Dataset"]))
    norm_data <- precisions_long %>% 
      inner_join(TPs_long) %>% 
      inner_join(positives_long)
    norm_data$`Normalisation method` <- factor(norm_data$`Normalisation method`,
                                               levels = c("Baseline",
                                                          "Intra-study"))
    
    # levels(norm_data$`Normalisation method`)
    refs_no <- c(length(refs[[1]][[1]]), length(refs[[2]][[1]]),
                 length(refs[[1]][[2]]), length(refs[[2]][[2]]))
    
    
    p1 <- norm_data %>% 
      ggplot(mapping = aes(x = `Normalisation method`,
                           y = Precision,
                           fill = Dataset)) +
      geom_col(position = "dodge",
               width = 0.4) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1),
            legend.direction = "horizontal") +
      annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
               x = 0.6, y = 1, hjust = 0)
    
    p2 <- norm_data %>%
      ggplot(mapping = aes(x = `Normalisation method`,
                           y = Positives,
                           fill = Dataset)) +
      geom_col(position = "dodge",
               width = 0.4) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1)) +
      annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
               x = 1.4, y = max(c(max(refs_no), max(norm_data$Positives))), hjust = 1) +
      ylab("Significant targets after swapping") +
      geom_segment(aes(x = c(0.5,0.5,1.5,1.5), xend = c(1.5,1.5,2.5,2.5),
                       y = c(refs_no[1],refs_no[2],refs_no[3],refs_no[4]),
                       yend = c(refs_no[1],refs_no[2],refs_no[3],refs_no[4])),
                   linetype = "dashed", color = c("#F8766D", "#00BFC4", "#F8766D", "#00BFC4"))
    
    legend <- get_legend(p1)
    g <- ggarrange(p1,p2,
                   ncol = 2,
                   common.legend = TRUE,
                   legend.grob = legend,
                   hjust = -2)
    title <- "Evaluation of normalisation method feasibility"
    annotate_figure(g, top = text_grob(title))
    # Noter til mig selv:
    # Sæt længde/bredde, så der kan være to af disse på et PP.
    # Nr. 2 skal være med FC_cutoff = 0.3.
    # Vi skal lige have en mere realistisk version af Limma (se One Note)
    # ^ Tjek. Det har vi nu. Hvordan ser det ud nu?
    filename <- paste0("figures/norm_perf_fc_", fc_cutoff)
    filename <- paste0(filename, "_1_5", ".png")
    ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 300, bg = "white")
}






################################################################################
# Designed to plot precision for regular setup vs. outliers removed
# Individual truths
for (plot in "letsgo"){
  # Define parameters
  fc_cutoff = 0.3
  
  cases <- list(
    list(
      log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # No normalisation (1-4)
      outrm_JYNR210602[outrm_meta_JYNR210602$Group == 3,] # Outlier removed
    ),
    list(
      log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
      log_JYNR220601[meta_mouse_JYNR220601$Group == 2,] # No normalisation (1-4)
    ))
  
  placebos <- list(
    list(
      log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
      log_JYNR220601[meta_mouse_JYNR220601$Group == 1,] # No normalisation (1-4)
    ),
    list(
      log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # No normalisation (1-4)
      outrm_JYNR210602[outrm_meta_JYNR210602$Group == 1,] # Outliers removed
    ))
  
  # Compute references (intra-study significant proteins)
  # Four references: Two for each study, one with intra-study normalisation (1-5) and one without (1-4)
  refs <- list(
    list(
      get_signif_targets(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                         group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
      get_signif_targets(soma1 = outrm_JYNR210602, meta1 = outrm_meta_JYNR210602,
                         group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff)),
    list(
      get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                         group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
      get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                         group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff)))
  
  
  precisions <- TPs <- positives <- matrix(nrow = 2, ncol = length(cases[[1]]))
  precisions <- TPs <- positives <- as.data.frame(precisions)
  rownames(precisions) <- rownames(TPs) <- rownames(positives) <- datasets <- c("JYNR210602", "JYNR220601")
  colnames(precisions) <- colnames(TPs) <- colnames(positives) <- c("Baseline",
                                                                    "Outliers removed")
  
  for (k in 1:2){ # for each dataset (JYNR210602 and JYNR220601)
    for (i in 1:dim(precisions)[2]){
      case <- cases[[k]][[i]]; placebo <- placebos[[k]][[i]]
      ref <- refs[[k]][[i]]
      positive <- get_signif_targets(soma1 = case, soma2 = placebo, fc_cutoff = fc_cutoff,
                                     meta_aligned = TRUE, sig_list_only = TRUE)
      TP <- sum(ref %in% positive)
      FP <- sum(!positive %in% ref)
      precision <- TP/(TP+FP)
      precisions[k,i] <- precision
      TPs[k,i] <- TP
      positives[k,i] <- length(positive)
    }
  }
  
  
  # Put into long format and join dfs
  precisions$Dataset <- TPs$Dataset <- positives$Dataset <- rownames(precisions)
  
  precisions_long <- precisions %>% 
    pivot_longer(names_to = "Normalisation method",
                 values_to = "Precision",
                 cols = colnames(precisions[,!colnames(precisions) == "Dataset"]))
  TPs_long <- TPs %>% 
    pivot_longer(names_to = "Normalisation method",
                 values_to = "TPs",
                 cols = colnames(positives[,!colnames(positives) == "Dataset"]))
  positives_long <- positives %>% 
    pivot_longer(names_to = "Normalisation method",
                 values_to = "Positives",
                 cols = colnames(positives[,!colnames(positives) == "Dataset"]))
  norm_data <- precisions_long %>% 
    inner_join(TPs_long) %>% 
    inner_join(positives_long)
  norm_data$`Normalisation method` <- factor(norm_data$`Normalisation method`,
                                             levels = c("Baseline",
                                                        "Outliers removed"))
  
  # levels(norm_data$`Normalisation method`)
  refs_no <- c(length(refs[[1]][[1]]), length(refs[[2]][[1]]),
               length(refs[[1]][[2]]), length(refs[[2]][[2]]))
  
  
  p1 <- norm_data %>% 
    ggplot(mapping = aes(x = `Normalisation method`,
                         y = Precision,
                         fill = Dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          legend.direction = "horizontal") +
    annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
             x = 0.6, y = 1, hjust = 0)
  
  p2 <- norm_data %>%
    ggplot(mapping = aes(x = `Normalisation method`,
                         y = Positives,
                         fill = Dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1)) +
    annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
             x = 1.4, y = max(c(max(refs_no), max(norm_data$Positives))), hjust = 1) +
    ylab("Significant targets after swapping") +
    geom_segment(aes(x = c(0.5,0.5,1.5,1.5), xend = c(1.5,1.5,2.5,2.5),
                     y = c(refs_no[1],refs_no[2],refs_no[3],refs_no[4]),
                     yend = c(refs_no[1],refs_no[2],refs_no[3],refs_no[4])),
                 linetype = "dashed", color = c("#F8766D", "#00BFC4", "#F8766D", "#00BFC4"))
  
  legend <- get_legend(p1)
  g <- ggarrange(p1,p2,
                 ncol = 2,
                 common.legend = TRUE,
                 legend.grob = legend,
                 hjust = -2)
  title <- "Evaluation of normalisation method feasibility"
  annotate_figure(g, top = text_grob(title))
  # Noter til mig selv:
  # Sæt længde/bredde, så der kan være to af disse på et PP.
  # Nr. 2 skal være med FC_cutoff = 0.3.
  # Vi skal lige have en mere realistisk version af Limma (se One Note)
  # ^ Tjek. Det har vi nu. Hvordan ser det ud nu?
  filename <- paste0("figures/norm_perf_fc_", fc_cutoff)
  filename <- paste0(filename, "_outrm", ".png")
  ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 300, bg = "white")
}





################################################################################
# Designed to plot precision for regular setup vs. outliers removed, version 2
# Individual truths
for (plot in "letsgo"){
  # Define parameters
  fc_cutoff = 0.3
  
  cases <- list(
    list(
      log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # No normalisation (1-4)
      outrm_JYNR210602[outrm_meta_JYNR210602$Group == 3,], # Outlier removed
      outrm2_JYNR210602[outrm_meta_JYNR210602$Group == 3,] # Outlier removed, pre-normalisation
    ),
    list(
      log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
      log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # No normalisation (1-4)
      outrm2_JYNR220601[outrm_meta_JYNR220601$Group == 2,] # Outlier removed, pre-normalisation
    ))
  
  placebos <- list(
    list(
      log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
      log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # No normalisation (1-4)
      outrm2_JYNR220601[outrm_meta_JYNR220601$Group == 1,] # Outlier removed, pre-normalisation
    ),
    list(
      log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # No normalisation (1-4)
      outrm_JYNR210602[outrm_meta_JYNR210602$Group == 1,], # Outliers removed
      outrm2_JYNR210602[outrm_meta_JYNR210602$Group == 1,] # Outliers removed
    ))
  
  # Compute references (intra-study significant proteins)
  # Four references: Two for each study, one with intra-study normalisation (1-5) and one without (1-4)
  refs <- list(
    list(
      get_signif_targets(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                         group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
      get_signif_targets(soma1 = outrm_JYNR210602, meta1 = outrm_meta_JYNR210602,
                         group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
      get_signif_targets(soma1 = outrm2_JYNR210602, meta1 = outrm_meta_JYNR210602,
                         group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff)),
    list(
      get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                         group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
      get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                         group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
      get_signif_targets(soma1 = outrm2_JYNR220601, meta1 = outrm_meta_JYNR220601,
                         group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff)))
  
  
  precisions <- TPs <- positives <- matrix(nrow = 2, ncol = length(cases[[1]]))
  precisions <- TPs <- positives <- as.data.frame(precisions)
  rownames(precisions) <- rownames(TPs) <- rownames(positives) <- datasets <- c("JYNR210602", "JYNR220601")
  colnames(precisions) <- colnames(TPs) <- colnames(positives) <- c("Baseline",
                                                                    "Outliers removed",
                                                                    "Outliers removed v2")
  
  for (k in 1:2){ # for each dataset (JYNR210602 and JYNR220601)
    for (i in 1:dim(precisions)[2]){
      case <- cases[[k]][[i]]; placebo <- placebos[[k]][[i]]
      ref <- refs[[k]][[i]]
      positive <- get_signif_targets(soma1 = case, soma2 = placebo, fc_cutoff = fc_cutoff,
                                     meta_aligned = TRUE, sig_list_only = TRUE)
      TP <- sum(ref %in% positive)
      FP <- sum(!positive %in% ref)
      precision <- TP/(TP+FP)
      precisions[k,i] <- precision
      TPs[k,i] <- TP
      positives[k,i] <- length(positive)
    }
  }
  
  
  # Put into long format and join dfs
  precisions$Dataset <- TPs$Dataset <- positives$Dataset <- rownames(precisions)
  
  precisions_long <- precisions %>% 
    pivot_longer(names_to = "Normalisation method",
                 values_to = "Precision",
                 cols = colnames(precisions[,!colnames(precisions) == "Dataset"]))
  TPs_long <- TPs %>% 
    pivot_longer(names_to = "Normalisation method",
                 values_to = "TPs",
                 cols = colnames(positives[,!colnames(positives) == "Dataset"]))
  positives_long <- positives %>% 
    pivot_longer(names_to = "Normalisation method",
                 values_to = "Positives",
                 cols = colnames(positives[,!colnames(positives) == "Dataset"]))
  norm_data <- precisions_long %>% 
    inner_join(TPs_long) %>% 
    inner_join(positives_long)
  norm_data$`Normalisation method` <- factor(norm_data$`Normalisation method`,
                                             levels = c("Baseline",
                                                        "Outliers removed",
                                                        "Outliers removed v2"))
  
  # levels(norm_data$`Normalisation method`)
  refs_no <- c(length(refs[[1]][[1]]), length(refs[[2]][[1]]),
               length(refs[[1]][[2]]), length(refs[[2]][[2]]),
               length(refs[[1]][[3]]), length(refs[[2]][[3]]))
  
  
  p1 <- norm_data %>% 
    ggplot(mapping = aes(x = `Normalisation method`,
                         y = Precision,
                         fill = Dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          legend.direction = "horizontal") +
    annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
             x = 0.6, y = 1, hjust = 0)
  
  p2 <- norm_data %>%
    ggplot(mapping = aes(x = `Normalisation method`,
                         y = Positives,
                         fill = Dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1)) +
    annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
             x = 1.4, y = max(c(max(refs_no), max(norm_data$Positives))), hjust = 1) +
    ylab("Significant targets after swapping") +
    geom_segment(aes(x = c(0.5,0.5,1.5,1.5,2.5,2.5), xend = c(1.5,1.5,2.5,2.5,3.5,3.5),
                     y = c(refs_no),
                     yend = c(refs_no)),
                 linetype = "dashed", color = rep(c("#F8766D", "#00BFC4"), length(cases[[1]])))
  
  legend <- get_legend(p1)
  g <- ggarrange(p1,p2,
                 ncol = 2,
                 common.legend = TRUE,
                 legend.grob = legend,
                 hjust = -2)
  title <- "Evaluation of normalisation method feasibility"
  annotate_figure(g, top = text_grob(title))
  # Noter til mig selv:
  # Sæt længde/bredde, så der kan være to af disse på et PP.
  # Nr. 2 skal være med FC_cutoff = 0.3.
  # Vi skal lige have en mere realistisk version af Limma (se One Note)
  # ^ Tjek. Det har vi nu. Hvordan ser det ud nu?
  filename <- paste0("figures/norm_perf_fc_", fc_cutoff)
  filename <- paste0(filename, "_outrm2", ".png")
  ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 300, bg = "white")
}




################################################################################
# Designed to plot precision for init, old, log, init_1_5, old_1_5, log_1_5
# Individual truths
for (plot in "letsgo"){
  # Define parameters
  fc_cutoff = 0.3
  
  norm_methods <- c("Initial",
                    "Two outliers removed",
                    "Outliers removed",
                    "MedNorm",
                    "Two outliers removed + medNorm",
                    "Outliers removed + medNorm")
  
  cases <- list(
    list(
      init_JYNR210602[init_meta_JYNR210602$Group == 3,], # Initial (outliers present)
      old_log_JYNR210602[old_meta_JYNR210602$Group == 3,], # Two outliers removed
      log_JYNR210602[meta_mouse_JYNR210602$Group == 3,], # Five outliers removed
      init_1_5_JYNR210602[init_meta_JYNR210602$Group == 3,], # Initial (outliers present)
      old_1_5_JYNR210602[old_meta_JYNR210602$Group == 3,], # Two outliers removed
      log_1_5_JYNR210602[meta_mouse_JYNR210602$Group == 3,] # Five outliers removed
    ),
    list(
      init_JYNR220601[init_meta_JYNR220601$Group == 2,], # Initial (outliers present)
      old_log_JYNR220601[old_meta_JYNR220601$Group == 2,], # Two outliers removed
      log_JYNR220601[meta_mouse_JYNR220601$Group == 2,], # Five outliers removed
      init_1_5_JYNR220601[init_meta_JYNR220601$Group == 2,], # Initial (outliers present)
      old_1_5_JYNR220601[old_meta_JYNR220601$Group == 2,], # Two outliers removed
      log_1_5_JYNR220601[meta_mouse_JYNR220601$Group == 2,] # Five outliers removed
    ))
  
  placebos <- list(
    list(
      init_JYNR220601[init_meta_JYNR220601$Group == 1,], # Initial (outliers present)
      old_log_JYNR220601[old_meta_JYNR220601$Group == 1,], # Two outliers removed
      log_JYNR220601[meta_mouse_JYNR220601$Group == 1,], # Five outliers removed
      init_1_5_JYNR220601[init_meta_JYNR220601$Group == 1,], # Initial (outliers present)
      old_1_5_JYNR220601[old_meta_JYNR220601$Group == 1,], # Two outliers removed
      log_1_5_JYNR220601[meta_mouse_JYNR220601$Group == 1,] # Five outliers removed
    ),
    list(
      init_JYNR210602[init_meta_JYNR210602$Group == 1,], # Initial (outliers present)
      old_log_JYNR210602[old_meta_JYNR210602$Group == 1,], # Two outliers removed
      log_JYNR210602[meta_mouse_JYNR210602$Group == 1,], # Five outliers removed
      init_1_5_JYNR210602[init_meta_JYNR210602$Group == 1,], # Initial (outliers present)
      old_1_5_JYNR210602[old_meta_JYNR210602$Group == 1,], # Two outliers removed
      log_1_5_JYNR210602[meta_mouse_JYNR210602$Group == 1,] # Five outliers removed
    ))
  
  # Compute references (intra-study significant proteins)
  # Four references: Two for each study, one with intra-study normalisation (1-5) and one without (1-4)
  refs <- list(
    list(
      get_signif_targets(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                         group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
      get_signif_targets(soma1 = outrm_JYNR210602, meta1 = outrm_meta_JYNR210602,
                         group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
      get_signif_targets(soma1 = outrm2_JYNR210602, meta1 = outrm_meta_JYNR210602,
                         group1 = 3, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff)),
    list(
      get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                         group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
      get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                         group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff),
      get_signif_targets(soma1 = outrm2_JYNR220601, meta1 = outrm_meta_JYNR220601,
                         group1 = 2, group2 = 1, sig_list_only = TRUE, fc_cutoff = fc_cutoff)))
  
  
  precisions <- TPs <- positives <- matrix(nrow = 2, ncol = length(cases[[1]]))
  precisions <- TPs <- positives <- as.data.frame(precisions)
  rownames(precisions) <- rownames(TPs) <- rownames(positives) <- datasets <- c("JYNR210602", "JYNR220601")
  colnames(precisions) <- colnames(TPs) <- colnames(positives) <- norm_methods
  
  for (k in 1:2){ # for each dataset (JYNR210602 and JYNR220601)
    for (i in 1:dim(precisions)[2]){
      case <- cases[[k]][[i]]; placebo <- placebos[[k]][[i]]
      if (k == 1){
        ref_placebo <- placebos[[2]][[i]]
      } else {
        ref_placebo <- placebos[[1]][[i]]
      }
      ref <- get_signif_targets(soma1 = case, soma2 = ref_placebo, fc_cutoff = fc_cutoff,
                                meta_aligned = TRUE, sig_list_only = TRUE)
      positive <- get_signif_targets(soma1 = case, soma2 = placebo, fc_cutoff = fc_cutoff,
                                     meta_aligned = TRUE, sig_list_only = TRUE)
      TP <- sum(ref %in% positive)
      FP <- sum(!positive %in% ref)
      precision <- TP/(TP+FP)
      precisions[k,i] <- precision
      TPs[k,i] <- TP
      positives[k,i] <- length(positive)
    }
  }
  
  
  # Put into long format and join dfs
  precisions$Dataset <- TPs$Dataset <- positives$Dataset <- rownames(precisions)
  
  precisions_long <- precisions %>% 
    pivot_longer(names_to = "Normalisation method",
                 values_to = "Precision",
                 cols = colnames(precisions[,!colnames(precisions) == "Dataset"]))
  TPs_long <- TPs %>% 
    pivot_longer(names_to = "Normalisation method",
                 values_to = "TPs",
                 cols = colnames(positives[,!colnames(positives) == "Dataset"]))
  positives_long <- positives %>% 
    pivot_longer(names_to = "Normalisation method",
                 values_to = "Positives",
                 cols = colnames(positives[,!colnames(positives) == "Dataset"]))
  norm_data <- precisions_long %>% 
    inner_join(TPs_long) %>% 
    inner_join(positives_long)
  norm_data$`Normalisation method` <- factor(norm_data$`Normalisation method`,
                                             levels = norm_methods)
  
  # levels(norm_data$`Normalisation method`)
  refs_no <- c(length(refs[[1]][[1]]), length(refs[[2]][[1]]),
               length(refs[[1]][[2]]), length(refs[[2]][[2]]),
               length(refs[[1]][[3]]), length(refs[[2]][[3]]))
  
  
  p1 <- norm_data %>% 
    ggplot(mapping = aes(x = `Normalisation method`,
                         y = Precision,
                         fill = Dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          legend.direction = "horizontal") +
    annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
             x = 0.6, y = 1, hjust = 0)
  
  p2 <- norm_data %>%
    ggplot(mapping = aes(x = `Normalisation method`,
                         y = Positives,
                         fill = Dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1)) +
    annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
             x = 1.4, y = max(c(max(refs_no), max(norm_data$Positives))), hjust = 1) +
    ylab("Significant targets after swapping") +
    geom_segment(aes(x = c(0.5,0.5,1.5,1.5,2.5,2.5), xend = c(1.5,1.5,2.5,2.5,3.5,3.5),
                     y = c(refs_no),
                     yend = c(refs_no)),
                 linetype = "dashed", color = rep(c("#F8766D", "#00BFC4"), length(cases[[1]])))
  
  legend <- get_legend(p1)
  g <- ggarrange(p1,p2,
                 ncol = 2,
                 common.legend = TRUE,
                 legend.grob = legend,
                 hjust = -2)
  title <- "Evaluation of normalisation method feasibility"
  annotate_figure(g, top = text_grob(title))
  # Noter til mig selv:
  # Sæt længde/bredde, så der kan være to af disse på et PP.
  # Nr. 2 skal være med FC_cutoff = 0.3.
  # Vi skal lige have en mere realistisk version af Limma (se One Note)
  # ^ Tjek. Det har vi nu. Hvordan ser det ud nu?
  filename <- paste0("figures/norm_perf_fc_", fc_cutoff)
  filename <- paste0(filename, "_outrm2", ".png")
  ggsave(filename, width = 4000, height = 2000, units = "px", dpi = 300, bg = "white")
}
