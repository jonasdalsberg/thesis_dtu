## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel)
setwd("~/JDJG/clean/")
source("scripts/SomaScan_Mouse/overlap_function.R")

## Load data and perform initial wrangling (this should be further automated for future use)
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")


## Bootstrap analysis
# Get bootstrap dataframes and stats (sig.targets, pre- and post-swap, precision, etc.)
for (run_bootstrap in "runthisline"){
  
  fc_cutoff <- 0.3
  
  bootstrap_dfs <- list(); bootstrap_stats <- list()
  for (dataset in c("JYNR210602", "JYNR220601")){
    print(dataset)
    # for (norm in c(FALSE)){
    # for (norm in c(FALSE, soma_cross, "limma")){
    for (norm in c(FALSE)){
      set.seed(1337)
      
      if (norm == "soma_cross"){
        soma1 <- "crossnorm_JYNR210602"; soma2 <- "crossnorm_JYNR220601"
        norm_text <- "Cross-normalised (SomaLogic)"
        print ("SomaLogic cross-norm")
      } else if (norm == "limma"){
        soma1 <- "limma_JYNR210602"; soma2 <- "limma_JYNR220601"
        norm_text <- "Batch effect correction (limma)"
        print ("limma")
      } else {
        soma1 <- "log_JYNR210602"; soma2 <- "log_JYNR220601"
        # soma1 <- "log_1_5_JYNR210602"; soma2 <- "log_1_5_JYNR220601"
        norm_text <- "No cross-normalisation"
        print("log")
      }
      soma1 <- get(paste0(soma1)); soma2 <- get(paste0(soma2))
      
      
      
      if (dataset == "JYNR210602"){
        case_1_orig <- soma1[meta_mouse_JYNR210602$Group == 3,]
        placebo_1_orig <- soma1[meta_mouse_JYNR210602$Group == 1,]
        case_2 <- soma2[meta_mouse_JYNR220601$Group == 2,]
        placebo_2 <- soma2[meta_mouse_JYNR220601$Group == 1,]
        case_n <- 7; placebo_n <- 5
        boot_n_c <- 15; boot_n_p <- 10
        boot_n_c <- 45; boot_n_p <- 30
      } else {
        case_1_orig <- soma2[meta_mouse_JYNR220601$Group == 2,]
        placebo_1_orig <- soma2[meta_mouse_JYNR220601$Group == 1,]
        case_2 <- soma1[meta_mouse_JYNR210602$Group == 3,]
        placebo_2 <- soma1[meta_mouse_JYNR210602$Group == 1,]
        case_n <- 15; placebo_n <- 10
        boot_n_c <- 7; boot_n_p <- 5
        boot_n_c <- 4; boot_n_p <- 3
      }
      
      it <- 100
      
      bootstraps <- list()
      bootstrap_data <- data.frame(sig.targets = rep(NA, it),
                                   sig.swap.1 = rep(NA, it),
                                   sig.swap.2 = rep(NA, it),
                                   TP.1 = rep(NA, it),
                                   TP.2 = rep(NA, it),
                                   precision.1 = rep(NA, it),
                                   precision.2 = rep(NA, it),
                                   fc_cutoff = rep(fc_cutoff, it))
      sig_tar_2 <- get_signif_targets(soma1 = case_2, soma2 = placebo_2,
                                      meta_aligned = TRUE, sig_list_only = TRUE,
                                      fc_cutoff = fc_cutoff)
      for (i in 1:it){
        print(i)
        case_i <- round(runif(boot_n_c, min = 1, max = case_n), 0)
        placebo_i <- round(runif(boot_n_p, min = 1, max = placebo_n), 0)
        # case_i <- 1:15; placebo_i <- 1:10
        case_1 <- case_1_orig[case_i,]; placebo_1 <- placebo_1_orig[placebo_i,]
        
        # Save bootstrapped dfs
        bootstrap <- list(case_1, placebo_1)
        bootstraps <- append(bootstraps, list(bootstrap))
        
        # Compute significant targets
        sig_tar_1 <- get_signif_targets(soma1 = case_1, soma2 = placebo_1,
                                        meta_aligned = TRUE, sig_list_only = TRUE,
                                        fc_cutoff = fc_cutoff)
        bootstrap_data$sig.targets[i] <- length(sig_tar_1)
        
        # Swap placebo groups and recompute significant targets
        sig_swap_1 <- get_signif_targets(soma1 = case_1, soma2 = placebo_2,
                                         meta_aligned = TRUE, sig_list_only = TRUE,
                                         fc_cutoff = fc_cutoff)
        sig_swap_2 <- get_signif_targets(soma1 = case_2, soma2 = placebo_1,
                                         meta_aligned = TRUE, sig_list_only = TRUE,
                                         fc_cutoff = fc_cutoff)
        bootstrap_data$sig.swap.1[i] <- length(sig_swap_1)
        bootstrap_data$sig.swap.2[i] <- length(sig_swap_2)
        
        # Compute precision
        bootstrap_data$TP.1[i] <- sum(sig_swap_1 %in% sig_tar_1)
        bootstrap_data$TP.2[i] <- sum(sig_swap_2 %in% sig_tar_2)
        bootstrap_data$precision.1[i] <- bootstrap_data$TP.1[i] / bootstrap_data$sig.swap.1[i]
        bootstrap_data$precision.2[i] <- bootstrap_data$TP.2[i] / bootstrap_data$sig.swap.2[i]
      }
      bootstrap_dfs <- append(bootstrap_dfs, list(bootstraps))
      bootstrap_stats <- append(bootstrap_stats, list(bootstrap_data))
      
      
      if (norm == "soma_cross"){
        crossnorm_bootstrap_dfs <- bootstrap_dfs
        crossnorm_bootstrap_stats <- bootstrap_stats
      } else if (norm == "limma"){
        limma_bootstrap_dfs <- bootstrap_dfs
        limma_bootstrap_stats <- bootstrap_stats
      } else {
        log_bootstrap_dfs <- bootstrap_dfs
        log_bootstrap_stats <- bootstrap_stats
      }
    }
  }
  save("bootstrap_n_test", "log_bootstrap_dfs", "log_bootstrap_stats",
       file = "~/JDJG/data/curated/Mouse_studies/bootstrap_temp.RData")
}

# # Bootstrapping complete!
# save("log_bootstrap_dfs", "log_bootstrap_stats",
#      file = "~/JDJG/data/curated/Mouse_studies/curated_bootstrap_mouse_1_5.RData")

# # Bootstrapping complete!
# save("log_bootstrap_dfs", "log_bootstrap_stats",
#      "crossnorm_bootstrap_dfs", "crossnorm_bootstrap_stats",
#      "limma_bootstrap_dfs", "limma_bootstrap_stats",
#      file = "~/JDJG/data/curated/Mouse_studies/curated_bootstrap_mouse.RData")
# save("log_bootstrap_dfs", "log_bootstrap_stats",
#      file = "~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_bootstrap_mouse.RData")







################################################################################
## Bootstrapping n test (figure out which Ns to use)
for (go in "go"){
  
  
  it <- 10000
  
  
  # fc_cutoff <- 0
  fc_cutoff <- 0.3
  
  # Define parameters
  n_tests_c <- c(7,15,30,45,90)
  n_tests_p <- c(5,10,20,30,60)
  # n_tests_c <- c(10000)
  # n_tests_p <- c(10000)
  
  skip_JYNR210602 <- FALSE
  # skip_JYNR210602 <- TRUE
  skip_JYNR220601 <- FALSE
  # skip_JYNR220601 <- TRUE
  
  # Don't touch the rest :-)
  l <- length(n_tests_c)
  m <- l * it * 2
  bootstrap_n_test <- data.frame(dataset = c(rep("JYNR210602", m/2), rep("JYNR220601", m/2)),
                                 n_c = rep(NA, m),
                                 n_p = rep(NA, m),
                                 it = rep(it,m),
                                 fc_cutoff = fc_cutoff,
                                 sig_tar = rep(NA, m),
                                 med_sig = rep(NA, m),avg_sig = rep(NA, m),
                                 sd = rep(NA, m),
                                 ci_low = rep(NA, m),
                                 ci_up = rep(NA, m),
                                 bootstrap = rep(NA,m))
  # Comment out:
  bootstrap_n_test <- bootstrap_n_test[1:(m-it),] # remove last empty n setting
  # from JYNR220601 (currently, 5 different ns are tested for JYNR210602 vs.
  # 4 different ns for JYNR220601)
  
  for (df_idx in c(0,m/2)){
    if (df_idx == 0){
      dataset <- "JYNR210602"
      case_df <- log_JYNR210602[meta_mouse_JYNR210602$Group == 3,]
      placebo_df <- log_JYNR210602[meta_mouse_JYNR210602$Group == 1,]
      case_n <- 7; placebo_n <- 5
    } else {
      dataset <- "JYNR220601"
      case_df <- log_JYNR220601[meta_mouse_JYNR220601$Group == 2,]
      placebo_df <- log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]
      case_n <- 15; placebo_n <- 10
      n_tests_c <- c(4,4,7,15)
      n_tests_p <- c(3,4,5,10)
    }
    
    for (i in 1:l){
      # Set seed
      set.seed(1337)
      
      
      # Skip the last n
      if (dataset == "JYNR220601" & i == l){
        break
      }
      if (dataset == "JYNR220601" & skip_JYNR220601 == TRUE){
        break
      }
      
      n_c <- n_tests_c[i]; n_p <- n_tests_p[i]
      idx <- (i-1)*it + df_idx
      for (j in 1:it){
        print(paste0(i, ", ", j, "."))
        # Ensure different samples are picked (not three identical samples for n = 3)
        case_i <- placebo_i <- c(1,1,1)
        while (sum(!case_i == case_i[1]) == 0 || sum(!placebo_i == placebo_i[1]) == 0){
          
          # If n is equal to the number of samples already present, simply keep the original samples
          case_i <- round(runif(n_c, min = 1, max = case_n), 0)
          placebo_i <- round(runif(n_p, min = 1, max = placebo_n), 0)
          #   if (n_c == case_n){
          #     case_i <- 1:case_n
          #   }
          #   if (n_p == placebo_n){
          #     placebo_i <- 1:placebo_n
          #   }
          #   
          #   # If n < group_n, sample randomly (should probably be without replacement)
          #   if (n_c < case_n){
          #     case_i <- round(runif(n, min = 1, max = case_n), 0)
          #   }
          #   if (n_p < placebo_n){
          #     placebo_i <- round(runif(n, min = 1, max = placebo_n), 0)
          #   }
          #   
          #   # If n > group_n, keep original samples and add until n == group_n
          #   if (n_c > case_n){
          #     n_new <- n - case_n
          #     case_i <- round(runif(n_new, min = 1, max = case_n), 0)
          #     case_i <- c(1:case_n, case_i)
          #   }
          #   if (n_p > placebo_n){
          #     n_new <- n - placebo_n
          #     placebo_i <- round(runif(n_new, min = 1, max = placebo_n), 0)
          #     placebo_i <- c(1:placebo_n, placebo_i)
          #   }
        }
        
        # Select samples
        case_temp <- case_df[case_i,]; placebo_temp <- placebo_df[placebo_i,]
        
        # Compute significant targets
        sig_targets <- get_signif_targets(soma1 = case_temp, soma2 = placebo_temp,
                                          meta_aligned = TRUE, sig_list_only = TRUE,
                                          fc_cutoff = fc_cutoff)
        
        # Store significant targets
        idx_j <- idx + j
        bootstrap_n_test$sig_tar[idx_j] <- length(sig_targets)
      }
      # Compute statistics for the current dataset and sample size (n)
      idx_1 <- idx + 1; int <- idx_1:idx_j
      bootstrap_n_test$n_c[int] <- n_c; bootstrap_n_test$n_p[int] <- n_p
      bootstrap_n_test$sd[int] <- sd(bootstrap_n_test$sig_tar[int])
      bootstrap_n_test$avg_sig[int] <- mean(bootstrap_n_test$sig_tar[int])
      bootstrap_n_test$med_sig[int] <- median(bootstrap_n_test$sig_tar[int]) 
      bootstrap_n_test$ci_low[int] <- bootstrap_n_test$avg_sig[idx_1] - 1.96 * bootstrap_n_test$sd[idx_1]
      bootstrap_n_test$ci_up[int] <- bootstrap_n_test$avg_sig[idx_1] + 1.96 * bootstrap_n_test$sd[idx_1]
      if (dataset == "JYNR210602"){
        if (n_c == 7 & n_p == 5){
          bootstrap_n_test$bootstrap[int] <- paste0("Original size, ", n_c, "/", n_p, "\nSomaMouse1")
        } else {
          bootstrap_n_test$bootstrap[int] <- paste0("Upsampled, ", n_c, "/", n_p)
        }
      } else {
        if (n_c == 15 & n_p == 10){
          bootstrap_n_test$bootstrap[int] <- paste0("Original size, ", n_c, "/", n_p, "\nSomaMouse2")
        } else {
          bootstrap_n_test$bootstrap[int] <- paste0("Downsampled, ", n_c, "/", n_p)
        }
      }
    }
  }
  bootstrap_n_test$bootstrap <- factor(bootstrap_n_test$bootstrap,
                                       levels = unique(bootstrap_n_test$bootstrap))
  save("bootstrap_n_test", "log_bootstrap_dfs", "log_bootstrap_stats",
       file = "~/JDJG/data/curated/Mouse_studies/bootstrap_temp.RData")
}

# bootstrap_n_test_fig7b_100$bootstrap <- NA
# for (i in 1:9){
#   k <- i*100; j <- k-99
#   int <- j:k
#   dataset <- unique(bootstrap_n_test_fig7b_100$dataset[int])
#   n_c <- unique(bootstrap_n_test_fig7b_100$n_c[int])
#   n_p <- unique(bootstrap_n_test_fig7b_100$n_p[int])
#   
#   
#   if (dataset == "JYNR210602"){
#     if (n_c == 7 & n_p == 5){
#       bootstrap_n_test_fig7b_100$bootstrap[int] <- paste0("Original size, ", n_c, "/", n_p, "\nSomaMouse1")
#     } else {
#       bootstrap_n_test_fig7b_100$bootstrap[int] <- paste0("Upsampled, ", n_c, "/", n_p)
#     }
#   } else {
#     if (n_c == 15 & n_p == 10){
#       bootstrap_n_test_fig7b_100$bootstrap[int] <- paste0("Original size, ", n_c, "/", n_p, "\nSomaMouse2")
#     } else {
#       bootstrap_n_test_fig7b_100$bootstrap[int] <- paste0("Downsampled, ", n_c, "/", n_p)
#     }
#   }
# }
# bootstrap_n_test_fig7b_100$bootstrap <- factor(bootstrap_n_test_fig7b_100$bootstrap,
#                                      levels = unique(bootstrap_n_test_fig7b_100$bootstrap))
# 
# unique(bootstrap_n_test_fig7b_100$n_c[801:900])
# bootstrap_n_test_fig7b_100$bootstrap[1:100]
# bootstrap_n_test_fig7b_100$testvar <- NA
# bootstrap_n_test_fig7b_100$testvar[1:100] <- paste0("Original size, ", n_c, "/", n_p, "\nSomaMouse1")
#<- bootstrap_n_test_fig7b_100$bootstrap[1:900]


# save("bootstrap_stats_fig7a", "bootstrap_dfs_fig7a",
#      "bootstrap_n_test_fig7b_100",
#      "bootstrap_stats_fig7c", "bootstrap_dfs_fig7c",
#      file = "~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/bootstrap_fig7.RData")

# save("bootstrap_n_test_10k",
#      file = "~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/bootstrap_10k.RData")


# JYNR210602_n_df_100 <- bootstrap_n_test
# bootstrap_100_5_15_21_28_35 <- plot_data
# bootstrap_100_8_15_22_30 <- plot_data
# bootstrap_100_up_down_realistic <- bootstrap_n_test
# bootstrap_levels <- c("Downsampled, 4/4",
#                       "Downsampled, 8/7",
#                       "Original size",
#                       "Upsampled, 15/10",
#                       "Upsampled, 20/20")
# bootstrap_100_up_down_realistic_2 <- bootstrap_n_test
# bootstrap_100_up_down_realistic_2$bootstrap <- NA
# bootstrap_100_up_down_realistic_2$bootstrap[c(1:it,(it*5+1):(it*6))] <- bootstrap_levels[3]
# bootstrap_100_up_down_realistic_2$bootstrap[c((it+1):(it*2))] <- bootstrap_levels[4]
# bootstrap_100_up_down_realistic_2$bootstrap[c((it*2+1):(it*3))] <- bootstrap_levels[5]
# bootstrap_100_up_down_realistic_2$bootstrap[c((it*3+1):(it*4))] <- bootstrap_levels[1]
# bootstrap_100_up_down_realistic_2$bootstrap[c((it*4+1):(it*5))] <- bootstrap_levels[2]
# bootstrap_100_up_down_realistic_2$bootstrap <- factor(bootstrap_100_up_down_realistic_2$bootstrap,
#                                                       levels = bootstrap_levels)

# save("bootstrap_100_5_15_21_28_35", "log_bootstrap_dfs", "log_bootstrap_stats",
#      file = "~/JDJG/data/curated/Mouse_studies/bootstrap_temp.RData")
# plot_data <- bootstrap_n_test
# plot_data$bootstrap[c(1:it,(it*(l+3)+1):(it*(l+4)))] <- "Original size"
# for (i in 2:l){
#   plot_data$bootstrap[c((it*(i-1)+1):(it*i))] <- paste0("Upsampled, ",
#                                                         plot_data$n_c[it*i],
#                                                         "/", plot_data$n_p[it*i])
# }
# for (i in (l+1):(l+3)){
#   plot_data$bootstrap[c((it*(i-1)+1):(it*i))] <- paste0("Downsampled, ",
#                                                         plot_data$n_c[it*i],
#                                                         "/", plot_data$n_p[it*i])
# }
# # plot_data[(it*5+1):(it*9),] <- bootstrap_5_up_down_realistic_2[(it*5+1):(it*9),]
# plot_data %>% 
#   ggplot(mapping = aes(x = dataset,
#                        y = sig_tar,
#                        fill = bootstrap)) +
#   geom_boxplot()+
#   ylab("Significant targets (intra-study)") +
#   xlab("Dataset") +
#   ggtitle(paste0("Bootstrapping test for different sampling sizes")) +
#   annotate("text", label = paste0("it = ", unique(plot_data$it), "\n",
#                                   "fc_cutoff = ", unique(plot_data$fc_cutoff)),
#            x = 0.5, y = max(c(plot_data$ci_up, plot_data$sig_tar)), hjust=0, vjust=1) +
#   scale_fill_manual(values = c("#00AFBB", "#00AFBB", "#00AFBB", "grey", "#bb0c00", "#bb0c00", "#bb0c00", "#bb0c00"))
# ggsave("figures/bootstrap_n_test_1000_fc_0_realistic_scale_boxplot.png", width = 2000, height = 1000, units = "px", dpi = 200, bg = "white")
# 
# 
# 
# # Plot the data
# plot_data <- bootstrap_100_up_down_realistic_2
# # plot_data <- bootstrap_n_test[bootstrap_n_test$dataset == "JYNR220601",]
# plot_data %>%
#   ggplot(mapping = aes(x = n_c,
#                        y = sig_tar,
#                        color = dataset)) +
#   geom_point() +
#   geom_line(aes(y = med_sig),
#             linetype = 2) +
#   geom_line(aes(y = avg_sig),
#             linetype = 1) +
#   geom_errorbar(mapping = aes(ymin = ci_low,
#                               ymax = ci_up,
#                               width = 0.2)) +
#   ggtitle(paste0("Bootstrapping test for different numbers of samples")) +
#   ylab("Significant targets (intra-study)") +
#   xlab("Number of samples per group (case vs. control), n") +
#   annotate("text", label = paste0("it = ", unique(plot_data$it)),
#            x = min(plot_data$n_c), y = max(c(plot_data$ci_up, plot_data$sig_tar)), hjust=0, vjust=1)
# # ggsave("figures/bootstrap_n_test_JYNR220601_100_4_8_15_30.png", width = 4000, height = 1200, units = "px", dpi = 200, bg = "white")
# 
# 
# plot_data %>% 
#   ggplot(mapping = aes(x = dataset,
#                        y = sig_tar,
#                        fill = bootstrap)) +
#   geom_boxplot()+
#   ylab("Significant targets (intra-study)") +
#   xlab("Dataset") +
#   ggtitle(paste0("Bootstrapping test for different numbers of samples")) +
#   annotate("text", label = paste0("it = ", unique(plot_data$it)),
#            x = 0.5, y = max(c(plot_data$ci_up, plot_data$sig_tar)), hjust=0, vjust=1) +
#   scale_fill_manual(values = c("#00AFBB", "#00AFBB", "grey", "#bb0c00", "#bb0c00"))
# ggsave("figures/bootstrap_n_test_100_realistic_scale_boxplot_2.png", width = 2000, height = 1000, units = "px", dpi = 200, bg = "white")






## FOR TOP100 BOOTSTRAPPING:
# Get bootstrap dataframes and stats (sig.targets, pre- and post-swap, precision, etc.)
for (run_bootstrap in "runthisline"){

  fc_cutoff <- 0.3
  it <- 20

  bootstrap_data <- data.frame(overlaps = rep(NA, it*2),
                               dataset = c(rep("SomaMouse1", it), rep("SomaMouse2", it)),
                               iteration = 1:(it*2))

  top100_JYNR210602 <- top100_JYNR220601 <- list()

  case_1_orig <- log_JYNR210602[meta_mouse_JYNR210602$Group == 3,]
  placebo_1_orig <- log_JYNR210602[meta_mouse_JYNR210602$Group == 1,]

  case_2_orig <- log_JYNR220601[meta_mouse_JYNR220601$Group == 2,]
  placebo_2_orig <- log_JYNR220601[meta_mouse_JYNR220601$Group == 1,]

  JYNR210602_sig <- get_signif_targets(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602,
                                       group1 = 3, group2 = 1)
  JYNR210602_sig$log2fc <- abs(JYNR210602_sig$log2fc)
  JYNR210602_sig <- JYNR210602_sig[order(JYNR210602_sig$log2fc, decreasing = TRUE),]
  top100_JYNR210602_orig <- JYNR210602_sig$seqid[1:100]

  JYNR220601_sig <- get_signif_targets(soma1 = log_JYNR220601, meta1 = meta_mouse_JYNR220601,
                                       group1 = 2, group2 = 1)
  JYNR220601_sig$log2fc <- abs(JYNR220601_sig$log2fc)
  JYNR220601_sig <- JYNR220601_sig[order(JYNR220601_sig$log2fc, decreasing = TRUE),]
  top100_JYNR220601_orig <- JYNR220601_sig$seqid[1:100]

  
  # JYNR210602
  set.seed(1337)
  for (i in 1:it){
    print(paste0("SomaMouse1, ", i))
    case_i <- round(runif(7, min = 1, max = 7), 0)
    placebo_i <- round(runif(5, min = 1, max = 5), 0)
    case_1 <- case_1_orig[case_i,]; placebo_1 <- placebo_1_orig[placebo_i,]
    sig_tar_1 <- get_signif_targets(soma1 = case_1, soma2 = placebo_1,
                                    meta_aligned = TRUE, fc_cutoff = fc_cutoff)
    sig_tar_1$log2fc <- abs(sig_tar_1$log2fc)
    sig_tar_1 <- sig_tar_1[order(sig_tar_1$log2fc, decreasing = TRUE),]
    top100_1 <- sig_tar_1$seqid[1:100]
    top100_JYNR210602 <- append(top100_JYNR210602,list(top100_1))
    overlap <- sum(top100_1 %in% top100_JYNR210602_orig)
    bootstrap_data$overlaps[i] <- overlap
  }
  
  
  # JYNR220601
  set.seed(1337)
  for (i in 1:it){
    print(paste0("SomaMouse2, ", i))
    case_i <- round(runif(15, min = 1, max = 15), 0)
    placebo_i <- round(runif(10, min = 1, max = 10), 0)
    case_2 <- case_2_orig[case_i,]; placebo_2 <- placebo_2_orig[placebo_i,]
    sig_tar_2 <- get_signif_targets(soma1 = case_2, soma2 = placebo_2,
                                    meta_aligned = TRUE, fc_cutoff = fc_cutoff)
    sig_tar_2$log2fc <- abs(sig_tar_2$log2fc)
    sig_tar_2 <- sig_tar_2[order(sig_tar_2$log2fc, decreasing = TRUE),]
    top100_2 <- sig_tar_2$seqid[1:100]
    top100_JYNR220601 <- append(top100_JYNR220601,list(top100_2))
    overlap <- sum(top100_2 %in% top100_JYNR220601_orig)
    bootstrap_data$overlaps[i+it] <- overlap
  }
  
}

# bootstrap_data_supfig16 <- bootstrap_data
# 
# save("bootstrap_data_supfig16",
#      file = "~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/bootstrap_supfig16.RData")

# Plot it

# sup_fig_14_data <- bootstrap_data_supfig14 <- bootstrap_data
# y1 <- median(sup_fig_14_data$overlaps[sup_fig_14_data$dataset == "SomaMouse1"])
# y2 <- median(sup_fig_14_data$overlaps[sup_fig_14_data$dataset == "SomaMouse2"])
# x1_start <- 0.5;
# x1_end <- dim(sup_fig_14_data)[1]/2 + 0.5;
# x2_start <- 0.5 + dim(sup_fig_14_data)[1]/2;
# x2_end <- dim(sup_fig_14_data)[1] + 0.5;
# sup_fig_14_data %>%
#   #filter(dataset == "JYNR210602") %>%
#   ggplot(mapping = aes(x = iteration,
#                        y = overlaps,
#                        fill = dataset)) +
#   geom_col(position = "dodge",
#            width = 0.4) +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 45,
#                                    hjust = 1),
#         plot.title = element_text(hjust = 0.5),
#         legend.direction = "vertical") +
#   annotate("text", label = paste0("Median = ", y1),
#            x = x1_start, y = 100, hjust = 0) +
#   annotate("text", label = paste0("Median = ", y2),
#            x = x2_end, y = 100, hjust = 1) +
#   scale_fill_manual(values = c("#FFA355", "#55B1FF")) +
#   geom_segment(x = x1_start, xend = x1_end, y = y1, yend = y1,
#                linetype = "dashed", color = "#FFA355") +
#   geom_segment(x = x2_start, xend = x2_end, y = y2, yend = y2,
#                linetype = "dashed", color = "#55B1FF") +
#   xlab("Bootstrap iteration") +
#   ylab("Overlaps between top100 genes") +
#   ggtitle("Top 100 genes by FC in original dataset versus boostrapped dataset")
# ggsave("figures/supfig14.png", width = 2000, height = 1000, units = "px", dpi = 200, bg = "white")


# bootstrap_data_supfig14 <- bootstrap_data
# Bootstrapping complete!
# save("bootstrap_data_supfig14",
#      file = "~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/bootstrap_supfig14.RData")









# Idea: Do a PCA and look for sample abnormalities.
# Is there a single sample outlier in JYNR220601?
# When I sample randomly, I do not see the big difference between the datasets (at least not as pronounced).