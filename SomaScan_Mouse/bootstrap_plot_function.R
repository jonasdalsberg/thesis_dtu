bootstrap_plot <- function(bootstrap_stats,
                           x_angle = 45, x_hjust = 1, x_vjust = 1,
                           short_y_lab = FALSE, n1 = c(7,5), n2 = c(15,10)){
  
  max_val <- max(unlist(bootstrap_stats))
  
  # Set theme
  theme_set(theme_classic() +
              theme(plot.title = element_text(hjust = 0.5, size = 12)))
  
  # Define standards
  stds <- c(479, 1924, 860, 1527)
  
  myboxplots <- list()
  for (i in 1:2){
    # Wrangle data
    bootstrap_data <- bootstrap_stats[[i]]
    it <- length(bootstrap_data$sig.targets)
    
    # change to long format
    data <- c("Intra-study 1",
              "Intra-study 2",
              "Swap 1",
              "Swap 2")
    df <- data.frame(`Differentially expressed targets` = c(rep(stds[1], it),
                                 rep(stds[2], it),
                                 bootstrap_data$sig.swap.1,
                                 bootstrap_data$sig.swap.2),
                     Data = factor(c(rep(data[1], it),
                                     rep(data[2], it),
                                     rep(data[3], it),
                                     rep(data[4], it)),
                                   levels = c(data)),
                     Dataset = c(rep("SomaMouse1", it),
                                 rep("SomaMouse2", it),
                                 rep("SomaMouse1", it),
                                 rep("SomaMouse2", it)))
    
    j <- 1
    if (i == 2){
      j <- j + 100
    }
    df$Differentially.expressed.targets[j:(j+99)] <- bootstrap_data$sig.targets
    
    if (short_y_lab == TRUE){
      y_lab <- "DE targets"
    } else {
      y_lab <- "Differentially expressed targets"
    }
    
    # Plot
    p <- df %>%
      ggplot(mapping = aes(x = Data,
                           y = Differentially.expressed.targets,
                           fill = Dataset)) +
      geom_boxplot() +
      ylim(c(0, max_val)) +
      scale_fill_manual(values = c("#FFA355", "#55B1FF")) +
      theme(axis.text.x = element_text(angle = x_angle,
                                       hjust = x_hjust, vjust = x_vjust)) +
      geom_segment(aes(x = 0.5, xend = 1.5,
                           y = stds[1], yend = stds[1]),
                   linetype = "dashed", color = "#FFA355") +
      geom_segment(aes(x = 1.5, xend = 2.5,
                       y = stds[2], yend = stds[2]),
                   linetype = "dashed", color = "#55B1FF") +
      geom_segment(aes(x = 2.5, xend = 3.5,
                       y = stds[3], yend = stds[3]),
                   linetype = "dashed", color = "#FFA355") +
      geom_segment(aes(x = 3.5, xend = 4.5,
                       y = stds[4], yend = stds[4]),
                   linetype = "dashed", color = "#55B1FF") +
      xlab("Data analysis") +
      ylab(y_lab)
    
    # Add iterations to the first plot
    if (i == 1){
      p <- p + 
        annotate("text", label = paste0("it = ", it),
                 x = 1.5, y = max_val, hjust=0, vjust=1)
    }
    
    myboxplots <- append(myboxplots, list(p))
  }
  
  # Create precision plots
  pre_stds <- c(0.423, 0.656)
  for (i in 1:2){
    # Wrangle data
    bootstrap_data <- bootstrap_stats[[i]]
    df <- data.frame(Precision = c(bootstrap_data$precision.1,
                                   bootstrap_data$precision.2),
                     Dataset = c(rep("SomaMouse1", it),
                                 rep("SomaMouse2", it)))
    
    p <- df %>% ggplot(mapping = aes(x = Dataset,
                                y = Precision,
                                fill = Dataset)) +
      geom_boxplot() +
      ylim(c(0, 1)) +
      scale_fill_manual(values = c("#FFA355", "#55B1FF")) +
      theme(axis.text.x = element_text(angle = x_angle,
                                       hjust = x_hjust, vjust = x_vjust)) +
      geom_segment(aes(x = 0.5, xend = 1.5,
                       y = pre_stds[1], yend = pre_stds[1]),
                   linetype = "dashed", color = "#FFA355") +
      geom_segment(aes(x = 1.5, xend = 2.5,
                       y = pre_stds[2], yend = pre_stds[2]),
                   linetype = "dashed", color = "#55B1FF")
    myboxplots <- append(myboxplots, list(p))
  }
  
  # Put the plots together
  legend <- get_legend(myboxplots[[1]])
  p1 <- ggarrange(myboxplots[[1]], myboxplots[[3]],
                  ncol = 2,
                  common.legend = TRUE,
                  legend = "none",
                  hjust = -2)
  p1 <- annotate_figure(p1, top = text_grob(paste0("Upsampling SomaMouse1, ", n1[1], "/", n1[2]), face = "bold"))

  p2 <- ggarrange(myboxplots[[2]], myboxplots[[4]],
                  ncol = 2,
                  common.legend = TRUE,
                  legend = "none",
                  hjust = -2)
  p2 <- annotate_figure(p2, top = text_grob(paste0("Downsampling SomaMouse2, ", n2[1], "/", n2[2]), face = "bold"))
  
  finalised_plots <- ggarrange(p1,p2,common.legend=TRUE,legend.grob=legend,legend="right",ncol = 2)
  
  return(finalised_plots)
}
