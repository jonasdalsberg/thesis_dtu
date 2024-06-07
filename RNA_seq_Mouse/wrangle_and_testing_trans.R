### Initialize document
# Load libraries
library(tidyverse); library(DESeq2)

# Set working directory
setwd("~/JDJG/clean/")

# Load functions
source("scripts/RNA_seq_Mouse/RNAseq_postprocessing_functions.R")
source("scripts/RNA_seq_Mouse/volcano_function.R")

# Load data
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/transcriptomics_mouse.RData")



################################################################################
#### Pre-process data
## Meta data
# TransMouse1
meta_trans1 <- gus053$meta_data[gus053$meta_data$group %in% c("Vehicle", "Chow vehicle"),]
meta_trans1$group <- as.character(meta_trans1$group)
meta_trans1[meta_trans1$group == "Vehicle",]$group <- "DIO-NASH Vehicle"
meta_trans1[meta_trans1$group == "Chow vehicle",]$group <- "Chow Vehicle"
meta_trans1$Dataset <- "TransMouse1"
meta_trans1$Treatment <- factor(meta_trans1$group,
                                levels = c("Chow Vehicle", 
                                           "DIO-NASH Vehicle"))
meta_trans1$Group <- meta_trans1$Model <- NA
meta_trans1$Group[meta_trans1$group == "DIO-NASH Vehicle"] <- 2;
meta_trans1$Group[meta_trans1$group == "Chow Vehicle"] <- 1;
meta_trans1$Model[meta_trans1$group == "DIO-NASH Vehicle"] <- "NASH"
meta_trans1$Model[meta_trans1$group == "Chow Vehicle"] <- "Healthy"
meta_trans1 <- meta_trans1[,-c(1,6)]
meta_trans1 <- meta_trans1[order(meta_trans1$Model),]
rownames(meta_trans1) <- meta_trans1$sample_id

# TransMouse2
meta_trans2 <- gus326$meta_data[gus326$meta_data$group %in% c("Chow Vehicle", "NASH Vehicle"),]
meta_trans2$group <- as.character(meta_trans2$group)
meta_trans2[meta_trans2$group == "NASH Vehicle",]$group <- "DIO-NASH Vehicle"
meta_trans2$Dataset <- "TransMouse2"
meta_trans2$Treatment <- factor(meta_trans2$group,
                                levels = c("Chow Vehicle", 
                                           "DIO-NASH Vehicle"))
meta_trans2$Group <- meta_trans2$Model <- NA
meta_trans2$Group[meta_trans2$group == "DIO-NASH Vehicle"] <- 2;
meta_trans2$Group[meta_trans2$group == "Chow Vehicle"] <- 1;
meta_trans2$Model[meta_trans2$group == "DIO-NASH Vehicle"] <- "NASH"
meta_trans2$Model[meta_trans2$group == "Chow Vehicle"] <- "Healthy"
meta_trans2 <- meta_trans2[,-c(1,6)]
meta_trans2 <- meta_trans2[order(meta_trans2$Model),]
rownames(meta_trans2) <- meta_trans2$sample_id

# TransMouse3
meta_trans3 <- gus257$meta_data[gus257$meta_data$group %in% c("Chow Vehicle", "DIO-NASH Vehicle"),]
meta_trans3$Dataset <- "TransMouse3"
meta_trans3$Treatment <- factor(meta_trans3$group,
                                levels = c("Chow Vehicle", 
                                           "DIO-NASH Vehicle"))
meta_trans3$Group <- meta_trans3$Model <- NA
meta_trans3$Group[meta_trans3$group == "DIO-NASH Vehicle"] <- 2;
meta_trans3$Group[meta_trans3$group == "Chow Vehicle"] <- 1;
meta_trans3$Model[meta_trans3$group == "DIO-NASH Vehicle"] <- "NASH"
meta_trans3$Model[meta_trans3$group == "Chow Vehicle"] <- "Healthy"
meta_trans3 <- meta_trans3[,-c(1,6)]
meta_trans3 <- meta_trans3[order(meta_trans3$Model),]
rownames(meta_trans3) <- meta_trans3$sample_id

## Expression data
# TransMouse1 (GUS053)
trans1 <- as.data.frame(t(gus053$exp_data$raw_counts))
colnames(trans1) <- unlist(trans1[1,])
anno_trans1 <- trans1[1:2,]
trans1 <- trans1[3:dim(trans1)[1],]
trans1 <- trans1[rownames(trans1) %in% meta_trans1$sample_id,] # keep only the two relevant groups
trans1 <- trans1[match(rownames(meta_trans1), rownames(trans1)),]
sum(!rownames(trans1) == rownames(meta_trans1)) # should be 0
trans1 <- mutate_all(trans1, function(x) as.numeric(as.character(x)))
# # Remove zero columns for either group
# rm <- c()
# for (i in 1:dim(trans1)[2]){
#   if (sum(trans1[meta_trans1$Group == 1,i]) == 0 || sum(trans1[meta_trans1$Group == 2,i]) == 0){
#     rm <- c(rm, i)
#   }
# }
# trans1 <- trans1[,-rm]

# TransMouse2 (GUS326)
trans2 <- as.data.frame(t(gus326$exp_data$raw_counts))
colnames(trans2) <- unlist(trans2[1,])
anno_trans2 <- trans2[1:2,]
trans2 <- trans2[3:dim(trans2)[1],]
trans2 <- trans2[rownames(trans2) %in% rownames(meta_trans2),] # keep only the two relevant groups
trans2 <- trans2[match(rownames(meta_trans2), rownames(trans2)),]
sum(!rownames(trans2) == rownames(meta_trans2)) # should be 0
trans2 <- mutate_all(trans2, function(x) as.numeric(as.character(x)))
# # Remove zero columns for either group
# rm <- c()
# for (i in 1:dim(trans2)[2]){
#   if (sum(trans2[meta_trans2$Group == 1,i]) == 0 || sum(trans2[meta_trans2$Group == 2,i]) == 0){
#     rm <- c(rm, i)
#   }
# }
# trans2 <- trans2[,-rm]

# TransMouse3 (GUS257)
trans3 <- as.data.frame(t(gus257$exp_data$raw_counts))
colnames(trans3) <- unlist(trans3[1,])
anno_trans3 <- trans3[1:2,]
trans3 <- trans3[3:dim(trans3)[1],]
trans3 <- trans3[rownames(trans3) %in% rownames(meta_trans3),] # keep only the two relevant groups
trans3 <- trans3[match(rownames(meta_trans3), rownames(trans3)),]
sum(!rownames(trans3) == rownames(meta_trans3)) # should be 0
trans3 <- mutate_all(trans3, function(x) as.numeric(as.character(x)))
# # Remove zero columns for either group
# rm <- c()
# for (i in 1:dim(trans3)[2]){
#   if (sum(trans3[meta_trans3$Group == 1,i]) == 0 || sum(trans3[meta_trans3$Group == 2,i]) == 0){
#     rm <- c(rm, i)
#   }
# }
# trans3 <- trans3[,-rm]

# # Align dfs
# trans1 <- trans1[,colnames(trans1) %in% colnames(trans2) & colnames(trans1) %in% colnames(trans3)]
# trans2 <- trans2[,colnames(trans2) %in% colnames(trans1)]
# trans3 <- trans3[,colnames(trans3) %in% colnames(trans1)]
# ens_id1 <- c(unlist(anno_trans1[1,]))
# ens_id2 <- c(unlist(anno_trans2[1,]))
# ens_id3 <- c(unlist(anno_trans3[1,]))
# not_in <- !ens_id2 %in% ens_id1
#### fix script (from here ^. incl. not_in to handle all datasets)
# ens_id <- c(ens_id1, ens_id2[not_in])
# gene1 <- c(unlist(anno_trans1[2,]))
# gene2 <- c(unlist(anno_trans2[2,]))
# gene <- c(gene1, gene2[not_in])
# anno_trans <- data.frame(ens_id = ens_id,
#                          gene = gene)
# anno_trans <- anno_trans[order(anno_trans[,1]),]
# rm("ens_id1", "ens_id2", "not_in", "ens_id",
#    "gene1", "gene2", "gene", "anno_trans1", "anno_trans2")
# 
# ## Already analysed data (dge)
# trans_dge1 <- gus053$dge_data$`Chow vehicle_vs_Vehicle`
# trans_dge2 <- gus326$dge_data$`Chow Vehicle_vs_NASH Vehicle`
# trans_dge3 <- gus257$dge_data$`DIO-NASH Vehicle_vs_Chow Vehicle`






# Plan:
# Skriv funktion, som kan tage to datasæt + to grupper og lave DGE
# Step 1:
# Lav wrangling ved brug af dine dfs, som er lig dette scripts wrangling
# Step 2:
# Lav wrangling, der kan splitte datasæt ad og samle dem igen med info om hvert sæt
# Step 3:
# Inkorporer det hele i en funktion
# Step 4:
# Udbyg funktion med følgende funktionaliteter:
# Swap TRUE/FALSE option
# Should add information into the output df that two different datasets have been compared
# And should change DESeqDataSetFromMatrix design from ~treatment to ~dataset
# Only necessary if control:control contrasts are being done


#### DE analysis: BATCH 1 #############################################
## Functions
# Auxilliary function
merge_dfs <- function(df1, df2){
  ex <- merge(t(df1), t(df2), by = "row.names", all = TRUE)
  rownames(ex) <- ex[,1]
  ex <- ex[-1]
  ex <- as.data.frame(t(ex))
  return (ex)
}

# Main function
dge_analysis <- function(df1, meta1, df2 = df1, meta2 = meta1, group1 = 2, group2 = 1,
                         meta_aligned = FALSE){
  df1 <- df1; df2 <- df2; meta1 <- meta1; meta2 <- meta2

  # Extract correct groups
  df1 <- df1[meta1$Group == group1,colnames(df1) %in% colnames(df2)]
  df2 <- df2[meta2$Group == group2,colnames(df2) %in% colnames(df1)]
  meta1 <- meta1[meta1$Group == group1,]
  meta2 <- meta2[meta2$Group == group2,]
  
  # Merge dataframes
  df <- merge_dfs(df1,df2)
  meta <- merge_dfs(meta1, meta2)
  rownames(meta) <- meta$sample_id
  
  # Get DESeq df
  if (group1 == group2){
    # If inter-study, use Dataset variable to contrast instead of Treatment
    dds <- DESeqDataSetFromMatrix(countData = round(t(df), 0),
                                  colData  = meta,
                                  design    = ~Dataset)
    grp1_name <- paste0(unique(meta1$Dataset))
    grp2_name <- paste0(unique(meta2$Dataset))
    contrast_list = list(c("Dataset", grp1_name, grp2_name))
    dds$Dataset <- relevel(dds$Dataset, ref = grp2_name)
    
  } else {
    dds <- DESeqDataSetFromMatrix(countData = round(t(df), 0),
                                  colData  = meta,
                                  design    = ~Treatment)
    
    grp1_name <- paste0(unique(meta1$Treatment))
    grp2_name <- paste0(unique(meta2$Treatment))
    contrast_list = list(c("Treatment", grp1_name, grp2_name))
    dds$Treatment <- relevel(dds$Treatment, ref = grp2_name)
  }
  
  # DESeq pipeline
  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)
  dds = nbinomWaldTest(dds)
  
  # Extract contrasts
  output_lists = extract_result_dfs(dds_model             = dds,
                                    extract_names         = contrast_list,
                                    intercept_or_contrast = "contrast",
                                    organism_name         = "mmusculus")
  output <- output_lists[[2]][[1]]
  
  return(output)
}



################################################################################
#### DGE analysis
### Analyse data
## Intra-study contrasts
dge_t1_intra <- dge_analysis(df1 = trans1, meta1 = meta_trans1)
dge_t2_intra <- dge_analysis(df1 = trans2, meta1 = meta_trans2)
dge_t3_intra <- dge_analysis(df1 = trans3, meta1 = meta_trans3)

## Inter-study (control:control)
# Chow vehicle
dge_t1_inter_chow_2 <- dge_analysis(df1 = trans1, meta1 = meta_trans1, df2 = trans2, meta2 = meta_trans2,
                                    group1 = 1, group2 = 1)
dge_t1_inter_chow_3 <- dge_analysis(df1 = trans1, meta1 = meta_trans1, df2 = trans3, meta2 = meta_trans3,
                              group1 = 1, group2 = 1)
dge_t2_inter_chow_1 <- dge_analysis(df1 = trans2, meta1 = meta_trans2, df2 = trans1, meta2 = meta_trans1,
                              group1 = 1, group2 = 1)
dge_t2_inter_chow_3 <- dge_analysis(df1 = trans2, meta1 = meta_trans2, df2 = trans3, meta2 = meta_trans3,
                              group1 = 1, group2 = 1)
dge_t3_inter_chow_1 <- dge_analysis(df1 = trans3, meta1 = meta_trans3, df2 = trans1, meta2 = meta_trans1,
                              group1 = 1, group2 = 1)
dge_t3_inter_chow_2 <- dge_analysis(df1 = trans3, meta1 = meta_trans3, df2 = trans2, meta2 = meta_trans2,
                              group1 = 1, group2 = 1)

# DIO-NASH vehicle
dge_t1_inter_nash_2 <- dge_analysis(df1 = trans1, meta1 = meta_trans1, df2 = trans2, meta2 = meta_trans2,
                                    group1 = 2, group2 = 2)
dge_t1_inter_nash_3 <- dge_analysis(df1 = trans1, meta1 = meta_trans1, df2 = trans3, meta2 = meta_trans3,
                                    group1 = 2, group2 = 2)
dge_t2_inter_nash_1 <- dge_analysis(df1 = trans2, meta1 = meta_trans2, df2 = trans1, meta2 = meta_trans1,
                                    group1 = 2, group2 = 2)
dge_t2_inter_nash_3 <- dge_analysis(df1 = trans2, meta1 = meta_trans2, df2 = trans3, meta2 = meta_trans3,
                                    group1 = 2, group2 = 2)
dge_t3_inter_nash_1 <- dge_analysis(df1 = trans3, meta1 = meta_trans3, df2 = trans1, meta2 = meta_trans1,
                                    group1 = 2, group2 = 2)
dge_t3_inter_nash_2 <- dge_analysis(df1 = trans3, meta1 = meta_trans3, df2 = trans2, meta2 = meta_trans2,
                                    group1 = 2, group2 = 2)

## Placebo swaps
dge_t1_swap_2 <- dge_analysis(df1 = trans1, meta1 = meta_trans1, df2 = trans2, meta2 = meta_trans2)
dge_t1_swap_3 <- dge_analysis(df1 = trans1, meta1 = meta_trans1, df2 = trans3, meta2 = meta_trans3)
dge_t2_swap_1 <- dge_analysis(df1 = trans2, meta1 = meta_trans2, df2 = trans1, meta2 = meta_trans1)
dge_t2_swap_3 <- dge_analysis(df1 = trans2, meta1 = meta_trans2, df2 = trans3, meta2 = meta_trans3)
dge_t3_swap_1 <- dge_analysis(df1 = trans3, meta1 = meta_trans3, df2 = trans1, meta2 = meta_trans1)
dge_t3_swap_2 <- dge_analysis(df1 = trans3, meta1 = meta_trans3, df2 = trans2, meta2 = meta_trans2)

## Compile into lists
dge_intra <- list(dge_t1_intra, dge_t2_intra, dge_t3_intra)
dge_inter_chow <- list(dge_t1_inter_chow_2, dge_t1_inter_chow_3, dge_t2_inter_chow_1,
                       dge_t2_inter_chow_3, dge_t3_inter_chow_1, dge_t3_inter_chow_2)
dge_inter_nash <- list(dge_t1_inter_nash_2, dge_t1_inter_nash_3, dge_t2_inter_nash_1,
                       dge_t2_inter_nash_3, dge_t3_inter_nash_1, dge_t3_inter_nash_2)
dge_swap <- list(dge_t1_swap_2, dge_t1_swap_3, dge_t2_swap_1,
                 dge_t2_swap_3, dge_t3_swap_1, dge_t3_swap_2)

## Save DGE dfs
save("dge_t1_intra", "dge_t2_intra", "dge_t3_intra",
     "dge_t1_inter_chow_2", "dge_t1_inter_chow_3", "dge_t2_inter_chow_1",
     "dge_t2_inter_chow_3", "dge_t3_inter_chow_1", "dge_t3_inter_chow_2",
     "dge_t1_inter_nash_2", "dge_t1_inter_nash_3", "dge_t2_inter_nash_1",
     "dge_t2_inter_nash_3", "dge_t3_inter_nash_1", "dge_t3_inter_nash_2",
     "dge_t1_swap_2", "dge_t1_swap_3", "dge_t2_swap_1",
     "dge_t2_swap_3", "dge_t3_swap_1", "dge_t3_swap_2",
     "dge_intra", "dge_inter_chow", "dge_inter_nash", "dge_swap",
     file = "~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/trans_dge.RData")

# ### Intra-study DE targets, placebo swap / precision, overlap
# ## Pre-process
# # Remove NAs (in the padj column)
# sig1_intra <- dge_t1_intra[!is.na(dge_t1_intra$padj),]
# sig2_intra <- dge_t2_intra[!is.na(dge_t2_intra$padj),]
# sig1_swap <- dge_t1_swap_2[!is.na(dge_t1_swap_2$padj),]
# sig2_swap <- dge_t2_swap_1[!is.na(dge_t2_swap_1$padj),]
# 
# # Filter by adjusted p-values
# sig1_intra <- sig1_intra$identifier[sig1_intra$padj < 0.05]
# sig2_intra <- sig2_intra$identifier[sig2_intra$padj < 0.05]
# sig1_swap <- sig1_swap$identifier[sig1_swap$padj < 0.05]
# sig2_swap <- sig2_swap$identifier[sig2_swap$padj < 0.05]
# 
# 
# ## Intra-study DGE and overlap btw. datasets
# length(sig1_intra) # 6748 DE genes (trans1)
# length(sig2_intra) # 6990 DE genes (trans2)
# sum(sig1_intra %in% sig2_intra) # 4750 overlapping
# sum(sig1_intra %in% sig2_intra) / length(sig1_intra) # 70.4% overlap (rel. to. trans1)
# sum(sig1_intra %in% sig2_intra) / length(sig2_intra) # 68.0% overlap (rel. to. trans2)
# 
# 
# # Placebo swapping and precision
# # TransMouse1
# length(sig1_swap) # 6336 DE genes (trans1, placebo-swapped)
# sum(sig1_swap %in% sig1_intra) # 3903 true positives (trans1)
# sum(sig1_swap %in% sig1_intra) / length(sig1_swap) # 57.8% precision (trans1)
# 
# # TransMouse2
# length(sig2_swap) # 9397 DE genes (trans2, placebo-swapped)
# sum(sig2_swap %in% sig2_intra) # 5765 true positives (trans2)
# sum(sig2_swap %in% sig2_intra) / length(sig2_swap) # 82.5% precision (trans2)
# 
# 
# ## Pre-process (adding in the third dataset)
# # Remove NAs (in the padj column)
# sig3_intra <- dge_t3_intra[!is.na(dge_t3_intra$padj),]
# sig1_swap_3 <- dge_t1_swap_3[!is.na(dge_t1_swap_3$padj),]
# sig2_swap_3 <- dge_t2_swap_3[!is.na(dge_t2_swap_3$padj),]
# sig3_swap_1 <- dge_t3_swap_1[!is.na(dge_t3_swap_1$padj),]
# sig3_swap_2 <- dge_t3_swap_2[!is.na(dge_t3_swap_2$padj),]
# 
# # Filter by adjusted p-values
# sig3_intra <- sig3_intra$identifier[sig3_intra$padj < 0.05]
# sig1_swap_3 <- sig1_swap_3$identifier[sig1_swap_3$padj < 0.05]
# sig2_swap_3 <- sig2_swap_3$identifier[sig2_swap_3$padj < 0.05]
# sig3_swap_1 <- sig3_swap_1$identifier[sig3_swap_1$padj < 0.05]
# sig3_swap_2 <- sig3_swap_2$identifier[sig3_swap_2$padj < 0.05]
# 
# 
# ## Intra-study DGE and overlap btw. datasets
# length(sig3_intra) # 6105 DE genes (trans3), similar to the other two datasets
# # Trans1 vs. trans3
# sum(sig1_intra %in% sig3_intra) # 4042 overlapping
# sum(sig1_intra %in% sig3_intra) / length(sig1_intra) # 59.9% overlap (rel. to. trans1)
# sum(sig1_intra %in% sig3_intra) / length(sig3_intra) # 66.2% overlap (rel. to. trans3)
# 
# # Trans2 vs. trans3
# sum(sig1_intra %in% sig3_intra) # 4293 overlapping
# sum(sig2_intra %in% sig3_intra) / length(sig2_intra) # 61.4% overlap (rel. to. trans2)
# sum(sig2_intra %in% sig3_intra) / length(sig3_intra) # 70.3% overlap (rel. to. trans3)
# 
# # All combined
# sum(sig1_intra %in% sig2_intra & sig1_intra %in% sig3_intra) # 3397
# sum(sig1_intra %in% sig2_intra & sig1_intra %in% sig3_intra) / length(sig1_intra) # 50.3% overlap (rel. to. trans1)
# sum(sig1_intra %in% sig2_intra & sig1_intra %in% sig3_intra) / length(sig2_intra) # 48.6% overlap (rel. to. trans2)
# sum(sig1_intra %in% sig2_intra & sig1_intra %in% sig3_intra) / length(sig3_intra) # 55.6% overlap (rel. to. trans3)
# 
# ## Placebo swapping and precision
# # Trans1 vs. trans3
# length(sig1_swap_3) # 8538 DE genes (trans1, placebo-swapped)
# sum(sig1_swap_3 %in% sig1_intra) # 4655 true positives (trans1)
# sum(sig1_swap_3 %in% sig1_intra) / length(sig1_swap_3) # 54.5% precision (trans1)
# 
# # Trans2 vs. trans3
# length(sig2_swap_3) # 9528 DE genes (trans2, placebo-swapped)
# sum(sig2_swap_3 %in% sig2_intra) # 5621 true positives (trans2)
# sum(sig2_swap_3 %in% sig2_intra) / length(sig2_swap_3) # 59.0% precision (trans2)
# 
# # Trans3 vs. trans1
# length(sig3_swap_1) # 9212 DE genes (trans3, placebo-swapped)
# sum(sig3_swap_1 %in% sig3_intra) # 4477 true positives (trans3)
# sum(sig3_swap_1 %in% sig3_intra) / length(sig3_swap_1) # 48.6% precision (trans3)
# 
# # Trans3 vs. trans2
# length(sig3_swap_2) # 7318 DE genes (trans3, placebo-swapped)
# sum(sig3_swap_2 %in% sig3_intra) # 3874 true positives (trans3)
# sum(sig3_swap_2 %in% sig3_intra) / length(sig3_swap_2) # 52.9% precision (trans3)
# 
# 
# ### Top100 overlap
# ## Pre-process
# # Sort by absolute FC
# top100_1 <- dge_t1_intra[order(abs(dge_t1_intra$log2FoldChange), decreasing = TRUE),]
# top100_2 <- dge_t2_intra[order(abs(dge_t2_intra$log2FoldChange), decreasing = TRUE),]
# top100_3 <- dge_t3_intra[order(abs(dge_t3_intra$log2FoldChange), decreasing = TRUE),]
# 
# # Remove NAs
# top100_1 <- top100_1[!is.na(top100_1$padj),]
# top100_2 <- top100_2[!is.na(top100_2$padj),]
# top100_3 <- top100_3[!is.na(top100_3$padj),]
# 
# # Filter by significance
# top100_1 <- top100_1[top100_1$padj < 0.05,]
# top100_2 <- top100_2[top100_2$padj < 0.05,]
# top100_3 <- top100_3[top100_3$padj < 0.05,]
# 
# # Filter by top100 and select ids
# top100_1 <- top100_1[1:100, colnames(top100_1) == "identifier"]
# top100_2 <- top100_2[1:100, colnames(top100_2) == "identifier"]
# top100_3 <- top100_3[1:100, colnames(top100_3) == "identifier"]
# 
# ## Compute top100 overlap
# sum(top100_1 %in% top100_2) # 44
# sum(top100_1 %in% top100_3) # 41
# sum(top100_2 %in% top100_3) # 37
# 
# 
# ### Volcano plots
# volcano_plot(dge_intra[[2]])
# volcano_plot(dge_intra[[2]])