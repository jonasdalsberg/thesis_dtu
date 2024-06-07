## Script to assess the effects of (cross-) normalization

## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel); library(scater); library(grid)

## Load data
# load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")

## Define function
# Bland-Altman plot function
bland_altman <- function(soma1 = log_JYNR210602, meta1 = meta_mouse_JYNR210602, case1 = 3, control1 = 1,
                         soma2 = log_JYNR220601, meta2 = meta_mouse_JYNR220601, case2 = 2, control2 = 1,
                         comparison = "case", normalisation = "", xlim = "", ylim = "",
                         datanames = c("JYNR210602", "JYNR220601"), title = "",
                         short_y_lab = FALSE, frontpage = FALSE){
  
  # Set theme
  theme_set(theme_classic() +
              theme(plot.title = element_text(hjust = 0.5, size = 12)))
  
  
  # Intended to produce a Bland-Altman plot for placebo:placebo or case:case comparisons
  # Use: Define soma1 and soma2, e.g. normal (log_JYNR210602), cross-normalized (norm_JYNR210602), etc.
  # Define if comparison is "case" or "control" ("placebo" also works)
  if (comparison == "case"){
    groups = c(case1, case2)
    comp_lab <- "case vs. case"
    plot_col <- "#BF3D35"
  } else if (comparison %in% c("control", "placebo")){
    groups = c(control1, control2)
    comp_lab <- "control vs. control"
    plot_col <- "#53BF2E"
  }
  
  if (frontpage == FALSE){
    line_col = "red"
  } else {
    line_col = "black"
    plot_col <- "#C51137"
  }
  
  # Split into case/control dfs
  df1 <- soma1[meta1$Group == groups[1],]
  df2 <- soma2[meta2$Group == groups[2],]
  
  # Compute median per aptamer
  med1 <- apply(df1, 2, median)
  med2 <- apply(df2, 2, median)
  
  # Compute means and differences
  df <- data.frame(soma1 = med1,
                   soma2 = med2)
  df$avg <- rowMeans(df)
  df$diff <- df$soma1 - df$soma2
  
  # Compute average difference and 95% CI
  mean_diff <- mean(df$diff)
  CI <- mean_diff + c(-1,1) * 1.96*sd(df$diff)
  
  # Define title
  if (title == ""){
    
    title <- paste0("Bland-Altman Plot, ", datanames[1], " vs. ", datanames[2],
                    ", ", comp_lab)
    if (!normalisation == ""){
      title <- paste0(title, ", ", normalisation)
    }
  }
  
  # Define y-axis label
  if (short_y_lab == TRUE){
    y_lab <- "Inter-study dif."
  } else {
    y_lab <- "Inter-study difference per aptamer"
  }
  
  # Plot it
  plt <- ggplot(df,
                mapping = aes(x = avg,
                              y = diff)) +
    geom_point(size=2, colour = plot_col, alpha = 0.5) +
    geom_hline(yintercept = mean_diff) +
    geom_hline(yintercept = CI[1],
               colour = line_col,
               linetype = "dashed") +
    geom_hline(yintercept = CI[2],
               colour = line_col,
               linetype = "dashed") +
    ggtitle(title) +
    ylab(y_lab) +
    xlab("Average cross-study RFU per aptamer")
  
  if (!xlim == ""){
    plt <- plt + xlim(xlim)
  }
  if (!ylim == ""){
    plt <- plt + ylim(ylim)
  }
  
  if (frontpage == TRUE){
    # Remove labels, numbers, titles, and change the colour scheme to bw
    plt <- plt +
      labs(x = "", y = "") +
      ggtitle("") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())
  }
  
  return(plt)
}

# Support function (xlim/ylim determination for Bland-Altman plotting)
determine_limits <- function(soma1_dfs = list(log_JYNR210602), meta1 = meta_mouse_JYNR210602, case1 = 3, control1 = 1,
                             soma2_dfs = list(log_JYNR220601), meta2 = meta_mouse_JYNR220601, case2 = 2, control2 = 1,
                             comp = "both"){
  # comp can be set to "case", "control", or "both"
  if (comp == "both"){
    comparisons <- list(c(case1, case2),
                        c(control1, control2))
  } else if (comp == "case"){
    comparisons <- list(c(case1, case2))
  } else if (comp == "control"){
    comparisons <- list(c(control1, control2))
  }
  
  x_mins <- x_maxs <- y_mins <- y_maxs <- list()
  for (i in 1:length(soma1_dfs)){
    for (groups in comparisons){
      
      # Get dfs
      soma1 <- soma1_dfs[[i]]
      soma2 <- soma2_dfs[[i]]
      
      # Split into case/control dfs
      df1 <- soma1[meta1$Group == groups[1],]
      df2 <- soma2[meta2$Group == groups[2],]
      
      # Compute median per aptamer
      med1 <- apply(df1, 2, median)
      med2 <- apply(df2, 2, median)
      
      # Compute means and differences
      df <- data.frame(soma1 = med1,
                       soma2 = med2)
      df$avg <- rowMeans(df)
      df$diff <- df$soma1 - df$soma2
      
      # Compute min/max
      x_mins <- append(x_mins, min(df$avg))
      x_maxs <- append(x_maxs, max(df$avg))
      y_mins <- append(y_mins, min(df$diff))
      y_maxs <- append(y_maxs, max(df$diff))
      
    }
  }
  limits <- c(min(unlist(x_mins)),
              max(unlist(x_maxs)),
              min(unlist(y_mins)),
              max(unlist(y_maxs)))
  return(limits)
}


# ################################################################################
# ## Plot it all
# # Designed to create Bland-Altman plots for case vs. case and control vs. control
# # for all three cross-normalization settings (including no normalization)
# for (all_in_one in "yes"){
#   norms <- c(FALSE, TRUE, "limma")
#   for (comp in c("case", "control")){
#     myplots <- list()
#     limits <- determine_limits(soma1_dfs = list(log_JYNR210602, norm_JYNR210602, limma_JYNR210602),
#                                soma2_dfs = list(log_JYNR220601, norm_JYNR220601, limma_JYNR220601),
#                                comp = comp)
#     xlim <- c(limits[1], limits[2]); ylim <- c(limits[3], limits[4])
#     for (norm in norms){
#       # Define dataframe
#       if (norm == TRUE){
#         soma1 <- "norm_JYNR210602"; soma2 <- "norm_JYNR220601"
#         norm <- "Cross-normalised (SomaLogic)"
#         print("cross-norm")
#       } else if (norm == "limma"){
#         soma1 <- "limma_JYNR210602"; soma2 <- "limma_JYNR220601"
#         norm <- "Batch effects corrected (limma)"
#         print ("limma")
#       } else {
#         soma1 <- "log_JYNR210602"; soma2 <- "log_JYNR220601"
#         norm <- "No normalisation"
#         print("log")
#       }
#       soma1 <- get(paste0(soma1)); soma2 <- get(paste0(soma2))
#       
#       # Plot
#       plt <- list(bland_altman(soma1 = soma1, soma2 = soma2, normalisation = norm,
#                           comparison = comp, title = norm,
#                           xlim = xlim, ylim = ylim))
#       myplots <- append(myplots, plt)
#     }
#     g <- ggarrange(plotlist = myplots,
#                    ncol = length(norms),
#                    nrow = 1,
#                    hjust = -2)
#     title <- paste0("Bland-Altman Plot, JYNR210602 vs. JYNR220601, ", comp, " vs. ", comp)
#     annotate_figure(g, top = text_grob(title))
#     filename <- paste0("~/JDJG/figures/bland_altman_", comp, ".png")
#     ggsave(filename, width = 4000, height = 1200, units = "px", dpi = 200, bg = "white")
#   }
# }
# 
# 
# 
# 
# 
# bland_altman(comparison = "control")
# ?element_text()
# 
# ggarrange(p1,p2)
# 
# 
# ggplot(case,
#        mapping = aes(x = avg,
#                      y = diff)) +
#   geom_point() +
#   xlim(c(0,20)) +
#   ggtitle("test") +
#   theme(plot.title = element_text(hjust = 0.5))
