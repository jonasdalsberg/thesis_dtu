



get_signif_targets <- function(soma1, meta1, soma2 = soma1, meta2 = meta1, group1, group2,
                           alpha = 0.05, mtc = "bh", log_t = TRUE, fc_cutoff = 0.3,
                           meta_aligned = FALSE, sig_list_only = FALSE){
  # For either two datasets or one datasets with two different treatment options,
  # this function computes log2(fc) and p-values for each comparison of SOMAmers
  # Option to use correction for multiple testing (Benjamini-Hochberg or Bonferroni)
  # Also provides a significance cutoff level for plotting.
  # Returns a dataframe.
  # Note: This function requires a SomaScanAnnotation file to be loaded.
  
  mtc <- stringr::str_to_lower(mtc)
  if (mtc == "bh"){
    mtc <- "BH"
  }
  
  if (meta_aligned == TRUE){
    # compute FC
    if (log_t == TRUE){
      fc <- sapply(1:ncol(soma1), function(x) mean(na.omit(as.vector(unlist(soma1[,x])))) - mean(na.omit(as.vector(unlist(soma2[,x])))))
    } else {
      fc <- log2(sapply(1:ncol(soma1), function(x) mean(na.omit(as.vector(unlist(soma1[,x])))) / mean(na.omit(as.vector(unlist(soma2[,x]))))))
    }
    
    # compute p-values (no correction for multiple testing)
    pval <- sapply(1:ncol(soma1), function(x) t.test(as.vector(unlist(soma2[,x])), as.vector(unlist(soma1[,x])) )$p.value)
    
  } else {
    # compute FC
    if (log_t == TRUE){
      fc <- sapply(1:ncol(soma1), function(x) mean(na.omit(as.vector(unlist(soma1[which(meta1$Group==group1),x])))) - mean(na.omit(as.vector(unlist(soma2[which(meta2$Group==group2),x])))))
    } else {
      fc <- log2(sapply(1:ncol(soma1), function(x) mean(na.omit(as.vector(unlist(soma1[which(meta1$Group==group1),x])))) / mean(na.omit(as.vector(unlist(soma2[which(meta2$Group==group2),x]))))))
    }
    
    # compute p-values (no correction for multiple testing)
    pval <- sapply(1:ncol(soma1), function(x) t.test(as.vector(unlist(soma2[which(meta2$Group == group2),x])), as.vector(unlist(soma1[which(meta1$Group == group1),x])) )$p.value)
  }
  
  # correct for multiple testing
  pval_adj <- p.adjust(pval, method = mtc)
  
  # Compute base mean
  base_mean_1 <- sapply(1:ncol(soma1), function(x) mean(as.vector(unlist(soma1[,x]))))
  base_mean_2 <- sapply(1:ncol(soma2), function(x) mean(as.vector(unlist(soma2[,x]))))
  
  # generate data frame
  volc <- data.frame(base_mean_1 = base_mean_1,
                     base_mean_2 = base_mean_2,
                     log2fc = fc,
                     pval = pval,
                     pval_adj = pval_adj)
  volc$difexpressed <- "NO"
  volc$difexpressed[volc$log2fc > fc_cutoff & volc$pval_adj < alpha] <- "UP"
  volc$difexpressed[volc$log2fc < -fc_cutoff & volc$pval_adj < alpha] <- "DOWN"
  volc$sign <- FALSE
  volc$sign[!volc$difexpressed == "NO"] <- TRUE
  volc$seqid = colnames(soma1); rownames(volc) <- colnames(soma1)
  
  # Check if SomaScan data or transcriptomics
  if (str_split(colnames(soma1)[1], "\\.")[[1]][1] == "seq"){
    volc$gene <- sapply(1:ncol(soma1), function(x) SomaScanAnnotation$Entrez.Gene.Name[SomaScanAnnotation$SeqId == colnames(soma1)[x]])
  } else {
    volc$gene <- sapply(1:ncol(soma1), function(x) anno_trans$gene[anno_trans$ens_id == colnames(soma1)[x]])
  }

    if (sig_list_only == TRUE){
    volc <- volc$seqid[volc$sign == TRUE]
  }
  
  return(volc)
}
