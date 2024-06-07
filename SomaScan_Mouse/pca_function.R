### LOAD LIBRARIES
# library(GEOquery)
library(limma)
library(tidyverse)
library(broom)
library(ggpubr)
library(ggfortify)
library(pheatmap)
library(RColorBrewer)
library(scales)

### DEFINE FUNCTIONS
# Boxplot function
plotbox <- function(df, x, y, legend_pos='none', col='none', xlab="", ylab=""){
  p <- df %>% 
    ggplot(mapping = aes_string(x = x,
                                y = y,
                                fill = x)) +
    geom_boxplot() +
    theme(legend.position=legend_pos) +
    labs(x = xlab, y = ylab) +
    theme_classic()
  if (!(col == 'none')){
    p <- p + scale_fill_brewer(palette=col)
  }
  return(p)
}

# Barplot function
plotbar <- function(df, x, col="none", title="none", xlab="", ylab=""){
  p <- df %>% 
    ggplot(mapping = aes_string(x = x,
                                fill = x)) +
    geom_bar(stat="count") +
    theme(legend.position="none") +
    labs(x = xlab, y = ylab) +
    theme_classic()
  if (!(col == "none")){
    p <- p + scale_fill_brewer(palette=col)
  }
  if (!(title == "none")){
    p <- p + ggtitle(title)
  }
  return(p)
}

pca_plot_1 <- function(df, meta, colour, shape, custom_cols, PC_axes = c(1,2),
                       title = " ", point_size = 5, short_labs = FALSE,
                       dataset){
  
  if (short_labs == TRUE){
    var_lab <- "% var)"
  } else {
    var_lab <- "% variance explained)"
  }
  
  var <- summary(df)$importance[2,]
  xlab <- paste0("PC", PC_axes[1], " (", round(var[PC_axes[1]]*100, 1), "% variance explained)")
  ylab <- paste0("PC", PC_axes[2], " (", round(var[PC_axes[2]]*100, 1), var_lab)
  
  p <- df %>% 
    augment(meta) %>% # add the meta data
    ggplot(mapping = aes_string(x = paste0(".fittedPC", PC_axes[1]),
                                y = paste0(".fittedPC", PC_axes[2]),
                                colour = colour,
                                shape = shape)) +
    geom_point(size = point_size) +
    theme_classic() +
    xlab(xlab) +
    ylab(ylab) +
    guides(colour = guide_legend(title=colour)) +
    custom_cols +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  if (length(dataset) == 1){
    if (dataset == "SomaMouse2"){
      
      p <- df %>% 
        augment(meta) %>% # add the meta data
        ggplot(mapping = aes_string(x = paste0(".fittedPC", PC_axes[1]),
                                    y = paste0(".fittedPC", PC_axes[2]),
                                    colour = colour,
                                    shape = shape)) +
        geom_point(size = point_size, shape = 17) +
        theme_classic() +
        xlab(xlab) +
        ylab(ylab) +
        guides(colour = guide_legend(title=colour)) +
        custom_cols +
        ggtitle(title) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    }
  }
  return(p)
}

# Plotting function for PCA
pcaplot <- function(data, colour = "Model", shape = "Dataset", groups = "all", title="none", legend="left", i=1,
                    center = TRUE, scale = TRUE, get_legend = FALSE, fig1 = FALSE, point_size = 5,
                    short_labs = FALSE, frontpage = FALSE){
  
  # Set theme
  theme_set(theme_classic() +
              theme(plot.title = element_text(hjust = 0.5)))
  
  df <- data[[1]]
  meta <- data[[2]]
  
  # Check if custom groups have been defined
  custom_groups <- FALSE
  if (!(length(groups) == 1)){
    custom_groups <- TRUE
  }
  if (length(groups) == 1){
    if (!(groups == "all")){
      custom_groups <- TRUE
    }
  }
  
  # Extract only relevant groups
  if (custom_groups == TRUE){
    include <- meta$Treatment %in% groups
    df <- df[include,]
    meta <- meta[include,]
  }
  
  # Compute principal components
  df <- prcomp(df, center = center, scale = scale)
  
  # Define specific colouring schemes for groups
  custom_cols <- scale_colour_discrete(hue_pal()(10))
  if (colour == "Model"){
    # print ("Colour == Model")
    custom_cols <- scale_colour_manual(values = c("#53BF2E", "#BF3D35"))
  }
  if (colour == "Treatment"){
    # print ("Colour == Treatment")
    if (length(unique(meta$Dataset)) == 2){
      # print ("Both datasets")
      custom_cols <- scale_colour_manual(values = c("#53BF2E", "#A9DF97",
                                                    "#BF3D35", "#CF6E68", "#DF9E9A", "#EFCFCD",
                                                    "#897E32"))
    } else {
      # print("Only one dataset")
      if (unique(meta$Dataset) == "SomaMouse1"){
        # print("JYNR210602")
        custom_cols <- scale_colour_manual(values = c("#53BF2E", "#A9DF97",
                                                      "#BF3D35", "#CF6E68", "#DF9E9A", "#EFCFCD"))
      } else if (unique(meta$Dataset) == "SomaMouse2"){
        # print("JYNR220601")
        custom_cols <- scale_colour_manual(values = c("#53BF2E",
                                                      "#BF3D35", "#897E32"))
      }
    }
  }
  if (!length(groups) == 1){
    # assume the vehicle groups have been chosen
    # (I should probably have made a dictionary for colour assignment instead)
    custom_cols <- scale_colour_manual(values = c("#53BF2E", "#BF3D35"))
  }
  
  
  # Create PCA plots
  if (frontpage == TRUE){
    # custom_cols <- scale_colour_manual(values = c("#53BF2E", "#A9DF97",
    #                                               "#BF3D35", "#CF6E68", "#DF9E9A", "#EFCFCD",
    #                                               "#897E32"))
    
    p1 <- pca_plot_1(df = df, meta = meta, colour = colour, shape = shape,
                     custom_cols = custom_cols, PC_axes = c(1,2),
                     point_size = point_size,
                     dataset = unique(meta$Dataset))
    # Remove labels, numbers, titles, and change the colour scheme
    p1 <- p1 +
      labs(colour = "", shape = "", x = "", y = "") +
      scale_colour_manual(values = c("#c51137", "#E2889B",
                                     "#001965", "#40538C", "#808CB2", "#BFC6D9",
                                     "#63154E")) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())
    return (p1)
  }
  
  p1 <- pca_plot_1(df = df, meta = meta, colour = colour, shape = shape,
                   custom_cols = custom_cols, PC_axes = c(1,2),
                   point_size = point_size, short_labs = short_labs,
                   dataset = unique(meta$Dataset))
  
  if (get_legend == TRUE){
    l <- get_legend(p1)
    return (l)
  }
  if (fig1 == TRUE){
    p2 <- pca_plot_1(df = df, meta = meta, colour = colour, shape = shape,
                     custom_cols = custom_cols, PC_axes = c(3,4), title = title,
                     point_size = point_size, short_labs = short_labs,
                     dataset = unique(meta$Dataset))
  } else {
    p2 <- pca_plot_1(df = df, meta = meta, colour = colour, shape = shape,
                     custom_cols = custom_cols, PC_axes = c(3,4),
                     point_size = point_size, short_labs = short_labs,
                     dataset = unique(meta$Dataset))
  }
  
  
  p <- ggarrange(p1,p2, common.legend = TRUE, legend=legend)
  # p1 <- pca_plot_1(df = df, meta = meta, colour = colour, custom_cols = custom_cols, PC_axes = c(1,2))
  # p2 <- pca_plot_1(df = df, meta = meta, colour = colour, custom_cols = custom_cols, PC_axes = c(2,3))
  # p3 <- pca_plot_1(df = df, meta = meta, colour = colour, custom_cols = custom_cols, PC_axes = c(3,4))
  # p4 <- pca_plot_1(df = df, meta = meta, colour = colour, custom_cols = custom_cols, PC_axes = c(4,5))
  # p <- ggarrange(p1,p2,p3,p4, common.legend = TRUE, legend=legend)
  if (!(title == "none") & !(fig1 == TRUE)){
    p <- annotate_figure(p, top = text_grob(title, face = "bold"))
  }
  return(p)
}

# Test plots
# pcaplot(log_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"),
#         title = "Mouse data, baseline, only overlapping control groups", legend = "left")
# pcaplot(log_data, groups = c("Chow Vehicle", "DIO-NASH Vehicle"),
#         title = "Mouse data, baseline, only overlapping control groups", legend = "left",
#         colour = "Treatment")
# pcaplot(log_data, title = "Mouse data, baseline, only overlapping control groups", legend = "left")
# pcaplot(log_data, title = "Mouse data, baseline, only overlapping control groups", legend = "left",
#         colour = "Treatment")
# pcaplot(list(log_JYNR210602, comb_meta[comb_meta$Dataset == "JYNR210602",]),
#         title = "JYNR210602", colour = "Treatment")
################################################################################

# Auxilliary function
merge_dfs <- function(df1, df2){
  ex <- merge(t(df1), t(df2), by = "row.names", all = TRUE)
  rownames(ex) <- ex[,1]
  ex <- ex[-1]
  ex <- as.data.frame(t(ex))
  return (ex)
}

# # Differential expression analysis function
# de_limma <- function(gse,meta,group){
#   ex <- exprs(gse)
#   sampleinfo <- meta %>% 
#     select(-dataset, -id)
#   design <- model.matrix(~0+sampleinfo[[group]])
#   names <- sort(unique(meta[[group]]))
#   colnames(design) <- c(names[1], names[2])
#   # Filtering lowly-expressed genes (cutoff = median)
#   cutoff <- median(ex)
#   # Subset genes expressed in more than 2 samples (maybe >2 is too few? 5 instead?)
#   is_expressed <- ex > cutoff
#   keep <- rowSums(is_expressed) > 2
#   table(keep) # check the ratio of remove:keep
#   ex <- ex[keep,]
#   # Coping with outliers (assigning weights to samples)
#   aw <- arrayWeights(ex,design)
#   # Fit a linear model to the data to estimate gene expression in group
#   fit1 <- lmFit(ex, design, weights = aw)
#   # Define the contrast for the differential analysis
#   contrasts <- paste0(colnames(design)[2], " - ", colnames(design)[1], ", levels=design")
#   contrasts <- makeContrastsFromString(contrasts, design = design)
#   fit2 <- contrasts.fit(fit1, contrasts)
#   # Apply empirical Bayes' step to get diff exp statistics and p-values
#   fit2 <- eBayes(fit2)
#   # Check how many differentially expressed genes we have
#   table(decideTests(fit2)) # 15937 (-1), 7809 (0), 12458 (1)
#   # Filter for significance
#   p_cutoff <- 0.05; fc_cutoff <- 1
#   merge_list <- topTable(fit2, number=Inf) %>%
#     tibble::rownames_to_column("ID") %>%
#     mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff)
#   return(merge_list)
# }

# # Auxilliary function
# makeContrastsFromString <- function(s, design){
#   eval(parse(text = paste("makeContrasts(", s, ")")))
# }