library(tidyverse); library(caret)
setwd("~/JDJG/clean")
source("scripts/SomaScan_Mouse/signif_targets_function.R")

get_sign_lists <- function(soma1, meta1, case1, control1,
                           soma2, meta2, control2,
                           alpha = 0.05, mtc = "none", log_t = TRUE, fc_cutoff = 0.3){
  # Provides a list of dataframes with significant targets for
    # Intra-study (case vs. control in soma 1)
    # Inter-study (control vs. control in soma 1 vs. soma 2)
    # Inter-study (case vs. control, placebo-swapped, i.e. case 1 vs. control 2)
  
  sign_lists <- list()
  for (scenario in list(list(soma1, meta1, case1, soma1, meta1, control1), # intra-study (case vs. control)
                        list(soma1, meta1, control1, soma2, meta2, control2), # inter-study (control vs. control)
                        list(soma1, meta1, case1, soma2, meta2, control2))){ # inter-study / placebo swab (case vs. control)
    df1 <- scenario[[1]]; meta_1 <- scenario[[2]]; group1 <- scenario[[3]]
    df2 <- scenario[[4]]; meta_2 <- scenario[[5]]; group2 <- scenario[[6]]
    sign_list <- get_signif_targets(soma1 = df1, meta1 = meta_1, group1 = group1,
                                    soma2 = df2, meta2 = meta_2, group2 = group2,
                                    alpha = alpha, mtc = mtc, log_t = log_t, fc_cutoff = fc_cutoff)

    sign_lists <- append(sign_lists, list(sign_list))
  }
  
  return (sign_lists)
}

compute_accuracy <- function(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602, case1 = 3, control1 = 1,
                              soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601, control2 = 1,
                              alpha = 0.05, mtc = "BH", log_t = TRUE, fc_cutoff = 0.3){
  
  # Get significant proteins
  sign_lists <- get_sign_lists(soma1, meta1, case1, control1,
                               soma2, meta2, control2,
                               alpha = alpha, mtc = mtc, log_t = log_t, fc_cutoff = fc_cutoff)
  
  
  # Confusion matrix test
  sign.prots <- factor(sign_lists[[1]]$sign, levels = c(TRUE, FALSE))
  sign.prots.swap <- factor(sign_lists[[3]]$sign, levels = c(TRUE, FALSE))
  cm <- confusionMatrix(sign.prots.swap, sign.prots, dnn = c("Original", "Placebo swapped"))
  metrics <- c(cm$byClass[[5]], cm$table[1], cm$table[3], cm$table[2], cm$byClass[[6]], 1-cm$byClass[[6]])
  # Precision, TP, FP, FN, Recall, 1 - Recall
  return(metrics)
}

plot_overlaps <- function(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602, case1 = 3, control1 = 1,
                         soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601, control2 = 1,
                         ref_soma, ref_meta, ref_case, ref_contr = FALSE,
                         alpha = 0.05, mtc = "BH", log_t = TRUE, fc_cutoff = 0.3,
                         title = "Confusion matrix, original vs. placebo swapped",
                         meta_aligned = FALSE){
  
  if (meta_aligned == FALSE){
    print("meta_aligned == FALSE")
    # Get significant proteins
    sign_lists <- get_sign_lists(soma1, meta1, case1, control1,
                                 soma2, meta2, control2,
                                 alpha = alpha, mtc = mtc, log_t = log_t, fc_cutoff = fc_cutoff)
    sign.prots <- factor(sign_lists[[1]]$sign, levels = c(TRUE, FALSE))
    sign.prots.swap <- factor(sign_lists[[3]]$sign, levels = c(TRUE, FALSE))
  } else {
    print("meta_aligned == TRUE")
    # Get swapped
    sign.prots.swap <- get_signif_targets(soma1 = soma1, soma2 = soma2, meta_aligned = TRUE,
                                          alpha = alpha, mtc = mtc, log_t = log_t, fc_cutoff = fc_cutoff)
    sign.prots.swap <- factor(sign.prots.swap$sign, levels = c(TRUE, FALSE))
  }
  if (ref_contr == TRUE){
    # Get reference
    print("Using reference for confusion matrix.")
    sign.prots <- get_signif_targets(soma1 = ref_soma, meta1 = ref_meta, group1 = ref_case, group2 = ref_contr,
                                     alpha = alpha, mtc = mtc, log_t = log_t, fc_cutoff = fc_cutoff)
    sign.prots <- factor(sign.prots$sign, levels = c(TRUE, FALSE))
  }

  print("got this far")
  cm <- confusionMatrix(sign.prots, sign.prots.swap, dnn = c("Original", "Placebo swapped"))
  cm$table[4] <- NA # set non-significant targets in both cases to NA to better compare positives
  plt <- as.data.frame(cm$table)
  plt$Original <- factor(plt$Original, levels = rev(levels(plt$Original)))
  p <- ggplot(plt,
         mapping = aes(Original, Placebo.swapped, fill = Freq)) +
    geom_tile() + geom_text(aes(label=Freq)) +
    labs(x = "Original", y = "Placebo swapped") +
    scale_x_discrete(labels=c("Not significant", "Significant")) +
    scale_y_discrete(labels=c("Significant", "Not significant")) +
    ggtitle(label = title) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.direction = "vertical")
}
