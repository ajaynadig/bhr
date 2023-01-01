BHR_meta <- function(ss_list, trait_list, annotations, num_blocks, genomewide_correction, fixed_genes,  
                     overdispersion, all_models, num_null_conditions, slope_correction, gwc_exclusion,
                     custom_variant_variances, intercept, start_time){
  
  if(length(ss_list) > 1 & length(trait_list) > 1){
    print(paste0("Aggregating ",length(trait_list)," traits over ",length(ss_list)," frequency-function bins."))
  } else if(length(ss_list) > 1 & length(trait_list) == 1){
    print(paste0("Aggregating 1 trait over ",length(ss_list)," frequency-function bins."))
  } else if(length(ss_list) == 1 & length(trait_list) > 1){
    print(paste0("Aggregating ",length(trait_list)," traits over 1 frequency-function bin."))
  } else {
    print(paste0("Aggregating 1 trait over 1 frequency-function bin."))
  }
  
  h2_bin <- list()
  jk_h2_bin <- list()
  fixed_h2se_bin <- list()
  counter = 1
  for (ss_file_name in ss_list){ #Aggregate over sum
    total_h2 <- list()
    jk_h2 <- list()
    fixed_h2se <- list()
    
    if (custom_variant_variances == FALSE){
      ss_file_name$variant_variance = 2*ss_file_name$AF*(1-ss_file_name$AF)
    } 
    
    for (trait in trait_list){ #Aggregate over mean
      ss_use <- ss_file_name[ss_file_name$phenotype_key == trait,]
      ss_use = ss_use %>% group_by(gene) %>% summarise(betas = list(beta),
                                                                         variant_variances = list(variant_variance),
                                                                         chromosome = first(chromosome),
                                                                         gene_position = first(gene_position),
                                                                         N = first(N),
                                                                         phenotype_key = first(phenotype_key),
                                                                         num_variants = length(variant_variance))
      
      ss_use$w_t_beta = sapply(1:nrow(ss_use), function(x) sum(ss_use$variant_variances[[x]]*ss_use$betas[[x]]))
      ss_use$burden_score = sapply(1:nrow(ss_use), function(x) sum(ss_use$variant_variances[[x]]))
      ss_use$cumulative_variance = sapply(1:nrow(ss_use), function(x) sum(ss_use$variant_variances[[x]]))
      ss_use$overdispersion = sapply(1:nrow(ss_use), function(x) sum(ss_use$variant_variances[[x]]^2)/sum(ss_use$variant_variances[[x]]))
      ss_use = ss_use[ss_use$burden_score > 0 & is.finite(ss_use$burden_score) & is.finite(ss_use$w_t_beta) & is.finite(ss_use$overdispersion),]

      object <- BHR_h2(ss_use, annotations, num_blocks, genomewide_correction, 
                       fixed_genes, output_jackknife_h2 = TRUE, 
                       overdispersion, all_models = TRUE, num_null_conditions, slope_correction, 
                       gwc_exclusion, intercept, start_time)

      total_h2[[trait]] <- object$mixed_model$heritabilities[1,ncol(object$mixed_model$heritabilities)]
      jk_h2[[trait]] <- object$subthreshold_genes$jackknife_h2[ncol(object$mixed_model$heritabilities),]
      fixed_h2se[[trait]] <- object$mixed_model$heritabilities[2,ncol(object$mixed_model$heritabilities)] - object$subthreshold_genes$heritabilities[2,ncol(object$mixed_model$heritabilities)]
    }
    h2_bin[[counter]] <- mean(unlist(total_h2))
    jk_h2_bin[[counter]] <- rowMeans(data.frame(jk_h2))
    fixed_h2se_bin[[counter]] <- unlist(fixed_h2se)
    print(paste0("Completed frequency-function bin ", counter))
    
    counter = counter + 1
    
  }
  h2_mixed_final <- Reduce(`+`, h2_bin)
  
  jk_h2_final_se <- sqrt(((num_blocks - 1)/num_blocks)*sum((Reduce(`+`, jk_h2_bin) - mean(Reduce(`+`, jk_h2_bin)))^2))
  fixed_h2se_bin_final <- sqrt(sum(unlist(fixed_h2se_bin)^2)) / length(trait_list)
  total_se <- sqrt(sum(jk_h2_final_se^2 + fixed_h2se_bin_final^2))
  
  return(list(
    aggregated_mixed_model_h2 = h2_mixed_final,
    aggregated_mixed_model_h2se = total_se))
}
