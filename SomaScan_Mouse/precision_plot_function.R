## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel)
source("scripts/SomaScan_Mouse/overlap_function.R")
source("scripts/SomaScan_Mouse/annotation_compass_function.R")


# Define function
precision_plot <- function(soma1_list, meta1_list, soma2_list, meta2_list,
                           norm_methods, fc_cutoff = 0.3, just_numbers = FALSE,
                           short_y_lab = FALSE, fixed_ref = FALSE,
                           x_angle = 45, x_hjust = 1, x_vjust = 1,
                           meta_aligned = FALSE, cases, placebos, apt_qual = "all"){
  
  if (!apt_qual == "all"){
    if (apt_qual == "medium" || apt_qual == "med"){
      apt_filt <- c("Medium_biosignal", "Medium_dilution")
    } else {
      apt_filt <- c(str_to_title(apt_qual))
    }
    seqids <- SomaScanAnnotation$SeqId[SomaScanAnnotation$mouse_quality %in% apt_filt & SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)]
  } else {
    seqids <- SomaScanAnnotation$SeqId[SomaScanAnnotation$SeqId %in% colnames(log_JYNR210602)]
  }
  
  
  # Extract case/control dfs and compute refs
  if (meta_aligned == FALSE){
    if (fixed_ref == TRUE){
      soma1 <- soma1_list[[1]]; soma2 <- soma2_list[[1]]
      ref1 <- get_signif_targets(soma1 = soma1[,colnames(soma1) %in% seqids],
                                 meta1 = meta1_list[[1]], fc_cutoff = fc_cutoff,
                                 group1 = 3, group2 = 1, sig_list_only = TRUE)
      ref2 <- get_signif_targets(soma1 = soma2[,colnames(soma2) %in% seqids],
                                 meta1 = meta2_list[[1]], fc_cutoff = fc_cutoff,
                                 group1 = 2, group2 = 1, sig_list_only = TRUE)
    }
    cases <- placebos <- refs <- list(list(), list())
    for (i in 1:length(soma1_list)){
      soma1 <- soma1_list[[i]]; soma2 <- soma2_list[[i]]
      meta1 <- meta1_list[[i]]; meta2 <- meta2_list[[i]]
      cases[[1]] <- append(cases[[1]], list(soma1[meta1$Group == 3,]))
      cases[[2]] <- append(cases[[2]], list(soma2[meta2$Group == 2,]))
      placebos[[1]] <- append(placebos[[1]], list(soma2[meta2$Group == 1,]))
      placebos[[2]] <- append(placebos[[2]], list(soma1[meta1$Group == 1,]))
      if (fixed_ref == FALSE){
        ref1 <- get_signif_targets(soma1 = soma1[,colnames(soma1) %in% seqids], meta1 = meta1, fc_cutoff = fc_cutoff,
                                   group1 = 3, group2 = 1, sig_list_only = TRUE)
        ref2 <- get_signif_targets(soma1 = soma2[,colnames(soma2) %in% seqids], meta1 = meta2, fc_cutoff = fc_cutoff,
                                   group1 = 2, group2 = 1, sig_list_only = TRUE)
      }
      refs[[1]] <- append(refs[[1]], list(ref1)); refs[[2]] <- append(refs[[2]], list(ref2))
    } 
  } else {  
    # Compute references (fixed)
    case <- cases[[1]][[1]]; placebo <- placebos[[2]][[1]]
    case <- case[,colnames(case) %in% seqids]; placebo <- placebo[,colnames(placebo) %in% seqids]
    ref1 <- get_signif_targets(soma1 = case, soma2 = placebo,
                               meta_aligned = TRUE, sig_list_only = TRUE)
    case <- cases[[2]][[1]]; placebo <- placebos[[1]][[1]]
    case <- case[,colnames(case) %in% seqids]; placebo <- placebo[,colnames(placebo) %in% seqids]
    ref2 <- get_signif_targets(soma1 = case, soma2 = placebo,
                               meta_aligned = TRUE, sig_list_only = TRUE)
    refs <- list(list(), list())
    for (i in 1:length(cases[[1]])){
      refs[[1]] <- append(refs[[1]], list(ref1)); refs[[2]] <- append(refs[[2]], list(ref2))
      
    }
  }
  
  # Set up dfs for precision calculation
  precisions <- TPs <- ref_count <- positives <- as.data.frame(matrix(nrow = 2, ncol = length(cases[[1]])))
  rownames(precisions) <- rownames(TPs) <- rownames(ref_count) <- rownames(positives) <- datasets <- c("SomaMouse1", "SomaMouse2")
  colnames(precisions) <- colnames(TPs) <- colnames(ref_count) <- colnames(positives) <- norm_methods
  
  # Compute precision, true positives, and total positives
  for (k in 1:2){ # for each dataset (JYNR210602 and JYNR220601)
    for (i in 1:dim(precisions)[2]){
      case <- cases[[k]][[i]]; placebo <- placebos[[k]][[i]]
      # Filter for SOMAmer quality
      case <- case[,colnames(case) %in% seqids]; placebo <- placebo[,colnames(placebo) %in% seqids]
      ref <- refs[[k]][[i]]
      positive <- get_signif_targets(soma1 = case, soma2 = placebo, fc_cutoff = fc_cutoff,
                                     meta_aligned = TRUE, sig_list_only = TRUE)
      TP <- sum(ref %in% positive)
      FP <- sum(!positive %in% ref)
      precision <- TP/(TP+FP)
      precisions[k,i] <- precision
      TPs[k,i] <- TP
      ref_count[k,i] <- length(ref)
      positives[k,i] <- length(positive)
    }
  }
  
  # Put into long format and join dfs
  precisions$Dataset <- TPs$Dataset <- ref_count$Dataset <- positives$Dataset <- rownames(precisions)
  
  if (length(cases[[1]]) == 1){
    # special case
    precisions_long <- precisions
    colnames(precisions_long)[1] <- "Precision"
    precisions_long$`Normalisation method` <- rep(norm_methods,2)
    
    TPs_long <- TPs
    colnames(TPs_long)[1] <- "TPs"
    TPs_long$`Normalisation method` <- rep(norm_methods,2)
    
    positives_long <- positives
    colnames(positives_long)[1] <- "Positives"
    positives_long$`Normalisation method` <- rep(norm_methods,2)
    
    ref_count_long <- ref_count
    colnames(ref_count_long)[1] <- "Ref_count"
    ref_count_long$`Normalisation method` <- rep(norm_methods,2)
    
  } else {
    precisions_long <- precisions %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "Precision",
                   cols = colnames(precisions[,!colnames(precisions) == "Dataset"]))
    TPs_long <- TPs %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "TPs",
                   cols = colnames(TPs[,!colnames(TPs) == "Dataset"]))
    positives_long <- positives %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "Positives",
                   cols = colnames(positives[,!colnames(positives) == "Dataset"]))
    ref_count_long <- ref_count %>% 
      pivot_longer(names_to = "Normalisation method",
                   values_to = "Ref_count",
                   cols = colnames(ref_count[,!colnames(ref_count) == "Dataset"]))
  }
  norm_data <- precisions_long %>% 
    inner_join(TPs_long) %>% 
    inner_join(positives_long) %>% 
    inner_join(ref_count_long)
  norm_data$`Normalisation method` <- factor(norm_data$`Normalisation method`,
                                             levels = norm_methods)
  
  # Return only metrics if just_numbers is TRUE
  if (just_numbers == TRUE){
    return (norm_data)
  }
  
  # Compute positions for dashed lines
  x_start <- x_end <- y_start <- c()
  i_ <- -0.5; norms <- length(cases[[1]])
  for (i in 1:norms){
    i_ <- i_ + 1; j <- i + norms
    x_start <- c(x_start, rep(i_, 2))
    x_end <- c(x_end, rep(i_+1, 2))
    y_start <- c(y_start, ref_count_long$Ref_count[i], ref_count_long$Ref_count[j])
  }
  y_end <- y_start
  
  # Define y-axis label
  if (short_y_lab == TRUE){
    y_lab <- "DE targets"
  } else {
    y_lab <- "DE targets post-swap"
  }
  
  # Plot
  p1 <- norm_data %>% 
    ggplot(mapping = aes(x = `Normalisation method`,
                         y = Precision,
                         fill = Dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = x_angle,
                                     hjust = x_hjust, vjust = x_vjust),
          legend.direction = "horizontal") +
    ylim(c(0,1)) +
    scale_fill_manual(values = c("#FFA355", "#55B1FF")) #+
    # annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
    #          x = 0.6, y = 1, hjust = 0)
  
  p2 <- norm_data %>%
    ggplot(mapping = aes(x = `Normalisation method`,
                         y = Positives,
                         fill = Dataset)) +
    geom_col(position = "dodge",
             width = 0.4) +
    theme_light() +
    theme(axis.text.x = element_text(angle = x_angle,
                                     hjust = x_hjust, vjust = x_vjust)) +
    # annotate("text", label = paste0("FC_cutoff = ", fc_cutoff),
    #          x = 1.4, y = max(c(max(norm_data$Ref_count), max(norm_data$Positives))), hjust = 1) +#
    ylab(y_lab) +
    scale_fill_manual(values = c("#FFA355", "#55B1FF")) +
    geom_segment(aes(x = x_start, xend = x_end,
                     y = y_start,
                     yend = y_end),
                 linetype = "dashed", color = rep(c("#FFA355", "#55B1FF"), length(cases[[1]])))
  
  legend <- get_legend(p1)
  g <- ggarrange(p1,p2,
                 ncol = 2,
                 common.legend = TRUE,
                 legend.grob = legend,
                 hjust = -2)
  title <- "Evaluation of swapping feasibility"
  if (!apt_qual == "all"){
    if (apt_qual == "low"){
      title <- paste0(title, ", low quality aptamers")
    }
    if (apt_qual == "medium" || apt_qual == "med"){
      title <- paste0(title, ", medium quality aptamers")
    }
    if (apt_qual == "high"){
      title <- paste0(title, ", high quality aptamers")
    }

  }
  annotate_figure(g, top = text_grob(title, face = "bold"))
}

# # Test function
# precision_plot(soma1_list = list(init_JYNR210602, old_log_JYNR210602, log_JYNR210602),
#                        meta1_list = list(init_meta_JYNR210602, old_meta_JYNR210602, meta_mouse_JYNR210602),
#                        soma2_list = list(init_JYNR220601, old_log_JYNR220601, log_JYNR220601),
#                        meta2_list = list(init_meta_JYNR220601, old_meta_JYNR220601, meta_mouse_JYNR220601),
#                        norm_methods = c("No outliers removed",
#                                         "Two outliers removed",
#                                         "Five outliers removed"))
