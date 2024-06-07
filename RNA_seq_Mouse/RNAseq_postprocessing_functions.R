library(DESeq2); library(gProfileR); library(gprofiler2)

extract_result_dfs = function(dds_model, extract_names, intercept_or_contrast = "intercept", organism_name = "mmusculus"){
  
  result_df_lists_all = list()
  result_df_lists     = list()
  
  for (a_name in extract_names){
    
    # Extract DEG
    if (intercept_or_contrast == "intercept"){
      result_df = results(dds_model, name=a_name)
    }
    
    if (intercept_or_contrast == "contrast"){
      result_df = results(dds_model, contrast = c(a_name[1], a_name[2], a_name[3]))
      a_name    = paste(a_name[2], a_name[3], sep = '_vs_')
    }
    
    # Save all data
    result_df_all                         = result_df
    result_df_all                         = result_df_all[order(result_df_all$padj, decreasing = F),]
    anno_df                               = gconvert(rownames(result_df_all), organism = organism_name, filter_na = F, mthreshold = 1)
    rownames(anno_df)                     = anno_df$input
    result_df_all$name                    = anno_df[rownames(result_df_all), 'name']
    result_df_all$identifier              = rownames(result_df_all)
    result_df_lists_all[[a_name]] = data.frame(result_df_all, stringsAsFactors = F)
    
    # Save only significant data
    result_df_sub       = result_df[!is.na(result_df$padj),]
    result_df_sub       = result_df_sub[result_df_sub$padj < 0.05,]
    result_df_sub       = result_df_sub[order(result_df_sub$padj, decreasing = F),]
    result_df_sub       = subset(result_df_sub, lfcSE < 1)
    result_df_sub$name  = anno_df[rownames(result_df_sub), 'name']
    result_df_sub       = data.frame(result_df_sub, stringsAsFactors = F)
    
    # save in list
    result_df_sub$identifier         = rownames(result_df_sub)
    result_df_lists[[a_name]] = result_df_sub
  }
  
  return(list("Significant_results" = result_df_lists,
              "All_results"         = result_df_lists_all))
}
