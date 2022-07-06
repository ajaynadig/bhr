BHR_meta <- function(ss_list_trait1, trait_list, annotations, num_blocks, genomewide_correction, fixed_genes,  
                     overdispersion, all_models, num_null_conditions, slope_correction, gwc_exclusion){
  
  if(length(ss_list_trait1) > 1 & length(trait_list) > 1){
    print(paste0("Aggregating ",length(trait_list)," traits over ",length(ss_list_trait1)," frequency-function bins."))
  } else if(length(ss_list_trait1) > 1 & length(trait_list) == 1){
    print(paste0("Aggregating 1 trait over ",length(ss_list_trait1)," frequency-function bins."))
  } else if(length(ss_list_trait1) == 1 & length(trait_list) > 1){
    print(paste0("Aggregating ",length(trait_list)," traits over 1 frequency-function bin."))
  } else {
    print(paste0("Aggregating 1 trait over 1 frequency-function bin."))
  }
  
  h2_bin <- list()
  jk_h2_bin <- list()
  fixed_h2se_bin <- list()
  counter = 1
  for (ss_file_name in ss_list_trait1){ #sum
    
    total_h2 <- list()
    jk_h2 <- list()
    fixed_h2se <- list()
    for (trait in trait_list){ #mean
      #n_bhr_trait <- head(ss_file_name[ss_file_name$phenotype_key == trait,"N"],1) #special line for meta if doing slope correction across traits
      object <- BHR_h2(ss_file_name[ss_file_name$phenotype_key == trait,], list(annotations), num_blocks, genomewide_correction, 
                       fixed_genes, output_jackknife_h2 = TRUE, 
                       overdispersion, all_models = TRUE, num_null_conditions, slope_correction, 
                       gwc_exclusion)
      total_h2[[trait]] <- object$mixed_model$heritabilities[1,ncol(annotations)]
      jk_h2[[trait]] <- object$subthreshold_genes$jackknife_h2[ncol(annotations),]
      fixed_h2se[[trait]] <- object$mixed_model$heritabilities[2,ncol(annotations)] - object$subthreshold_genes$heritabilities[2,ncol(annotations)]
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
BHR_meta <- function(ss_list_trait1, trait_list, annotations, num_blocks, genomewide_correction, fixed_genes,  
                     overdispersion, all_models, num_null_conditions, slope_correction, gwc_exclusion){
  
  if(length(ss_list_trait1) > 1 & length(trait_list) > 1){
    print(paste0("Aggregating ",length(trait_list)," traits over ",length(ss_list_trait1)," frequency-function bins."))
  } else if(length(ss_list_trait1) > 1 & length(trait_list) == 1){
    print(paste0("Aggregating 1 trait over ",length(ss_list_trait1)," frequency-function bins."))
  } else if(length(ss_list_trait1) == 1 & length(trait_list) > 1){
    print(paste0("Aggregating ",length(trait_list)," traits over 1 frequency-function bin."))
  } else {
    print(paste0("Aggregating 1 trait over 1 frequency-function bin."))
  }
  
  h2_bin <- list()
  jk_h2_bin <- list()
  fixed_h2se_bin <- list()
  counter = 1
  for (ss_file_name in ss_list_trait1){ #sum
    
    total_h2 <- list()
    jk_h2 <- list()
    fixed_h2se <- list()
    for (trait in trait_list){ #mean
      #n_bhr_trait <- head(ss_file_name[ss_file_name$phenotype_key == trait,"N"],1) #special line for meta if doing slope correction across traits
      object <- BHR_h2(ss_file_name[ss_file_name$phenotype_key == trait,], list(annotations), num_blocks, genomewide_correction, 
                       fixed_genes, output_jackknife_h2 = TRUE, 
                       overdispersion, all_models = TRUE, num_null_conditions, slope_correction, 
                       gwc_exclusion)
      total_h2[[trait]] <- object$mixed_model$heritabilities[1,ncol(annotations)]
      jk_h2[[trait]] <- object$subthreshold_genes$jackknife_h2[ncol(annotations),]
      fixed_h2se[[trait]] <- object$mixed_model$heritabilities[2,ncol(annotations)] - object$subthreshold_genes$heritabilities[2,ncol(annotations)]
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
