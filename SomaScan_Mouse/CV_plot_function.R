## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel)

## CV function
cv_plot <- function(soma1, soma2, meta1, meta2, text = "",
                    only_cont = FALSE, only_case = FALSE, x_max = FALSE,
                    custom_title = FALSE, get_legend = FALSE, fig1 = FALSE,
                    hide_legend = FALSE, frontpage = FALSE){
  # Set theme
  theme_set(theme_classic() +
              theme(plot.title = element_text(hjust = 0.5)))
  
  
  cvs <- matrix(nrow = dim(soma1)[2],
                ncol = 4)
  cvs[,1] <- sapply(1:ncol(soma1), function(x) sd(soma1[meta1$Group == 1,x]) / mean(soma1[meta1$Group == 1,x]) * 100)
  cvs[,2] <- sapply(1:ncol(soma2), function(x) sd(soma2[meta2$Group == 1,x]) / mean(soma2[meta2$Group == 1,x]) * 100)
  cvs[,3] <- sapply(1:ncol(soma1), function(x) sd(soma1[meta1$Group == 3,x]) / mean(soma1[meta1$Group == 3,x]) * 100)
  cvs[,4] <- sapply(1:ncol(soma2), function(x) sd(soma2[meta2$Group == 2,x]) / mean(soma2[meta2$Group == 2,x]) * 100)
  cvs <- as.data.frame(cvs); names(cvs) <- c("JYNR210602_placebo", "JYNR220601_placebo", "JYNR210602_case", "JYNR220601_case")
  
  if (get_legend == TRUE){
    p1 <- cvs %>% ggplot() +
      geom_histogram(aes(x = JYNR210602_placebo, fill = "SomaMouse1"), alpha = 0.8) +
      geom_histogram(aes(x = JYNR220601_placebo, fill = "SomaMouse2"), alpha = 0.8) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.direction = "vertical") +
      scale_fill_manual(values = c("#FFA355", "#55B1FF"))
    l <- get_legend(p1)
    return (l)
  }
  
  
  
  xmin <- min(cvs); xmax <- max(cvs)
  if (!x_max == FALSE){
    xmax <- x_max
  }
  
  bins = 30
  if (frontpage == TRUE){
    bins = 30
  }
  
  if (custom_title == TRUE){
    placebo_title <- case_title <- text
  } else {
    if (!text == ""){
      title_text <- paste0(", ", text)
    } else {
      title_text <- ""
    }
    # placebo_title <- paste0("Placebo vs. placebo", title_text)
    placebo_title <- "Placebo vs. placebo"
    # case_title <- paste0("Case vs. case", title_text)
    case_title <- "Case vs. case"
  }
  
  if (fig1 == TRUE){
    placebo_title <- " "
  }
  
  p1 <- cvs %>% ggplot() +
    geom_histogram(aes(x = JYNR210602_placebo, fill = "SomaMouse1"), alpha = 0.8, bins = bins) +
    geom_histogram(aes(x = JYNR220601_placebo, fill = "SomaMouse2"), alpha = 0.8, bins = bins) +
    xlab("Coefficient of variation (%) per aptamer") +
    ggtitle(placebo_title) +
    xlim(c(xmin,xmax)) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.direction = "horizontal") +
    scale_fill_manual(values = c("#FFA355", "#55B1FF"))
  
  p2 <- cvs %>% ggplot() +
    geom_histogram(aes(x = JYNR210602_case, fill = "SomaMouse1"), alpha = 0.8, bins = bins) +
    geom_histogram(aes(x = JYNR220601_case, fill = "SomaMouse2"), alpha = 0.8, bins = bins) +
    xlab("Coefficient of variation (%) per aptamer") +
    ggtitle(case_title) +
    xlim(c(xmin,xmax)) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.direction = "horizontal") +
    scale_fill_manual(values = c("#FFA355", "#55B1FF"))
  
  if (frontpage == TRUE){
    # Remove labels, numbers, titles, and change the colour scheme
    p1 <- p1 +
      labs(x = "", y = "") +
      ggtitle("") +
      scale_fill_manual(values = c("#001965", "#c51137")) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()) +
      xlim(c(0,17))
    return(p1)
  }
  
  if (only_cont == TRUE){
    p1
  } else if (only_case == TRUE){
    p2
  } else {
    legend <- get_legend(p1)
    if (hide_legend == TRUE){
      g <- ggarrange(p1,p2,
                     legend = "none",
                     hjust = -2) 
    } else {
      g <- ggarrange(p1,p2,
                     common.legend = TRUE,
                     legend.grob = legend,
                     hjust = -2) 
    }
    annotate_figure(g, top = text_grob(text, face = "bold"))
    # return(g)
  }
}