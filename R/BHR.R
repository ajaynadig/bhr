BHR <- function(mode = NULL,
                trait1_sumstats = NULL,
                trait2_sumstats = NULL,
                annotations,
                num_blocks = 100,
                genomewide_correction = FALSE,
                fixed_genes = NULL,
                output_jackknife_h2 = FALSE,
                output_jackknife_rg = FALSE,
                ss_list_trait1 = NULL,
                ss_list_trait2 = NULL,
                trait_list = NULL,
                overdispersion = FALSE,
                all_models = FALSE,
                slope_correction = NULL,
                num_null_conditions = 5,
                gwc_exclusion = NULL,
                custom_weights = NULL,
                null_stats = FALSE,
                intercept = TRUE,
                custom_variant_variances = FALSE) {

  #Calculate variance variances from AF, if not provided
  if (custom_variant_variances == FALSE){
    if (is.null(trait2_sumstats)){
      trait1_sumstats$variant_variance = 2*trait1_sumstats$AF*(1-trait1_sumstats$AF)
    } else {
      trait1_sumstats$variant_variance = 2*trait1_sumstats$AF*(1-trait1_sumstats$AF)
      trait2_sumstats$variant_variance = 2*trait2_sumstats$AF*(1-trait2_sumstats$AF)
    }
  } 
  
    #Aggregate variant-row summary statistics into gene-row summary statistics
    if (is.null(trait2_sumstats)){
      trait1_sumstats = trait1_sumstats %>% group_by(gene) %>% summarise(betas = list(beta),
                                                                         variant_variances = list(variant_variance),
                                                                         chromosome = first(chromosome),
                                                                         gene_position = first(gene_position),
                                                                         N = first(N),
                                                                         phenotype_key = first(phenotype_key))
    } else{
      trait1_sumstats = trait1_sumstats %>% group_by(gene) %>% summarise(betas = list(beta),
                                                                         variant_variances = list(variant_variance),
                                                                         chromosome = first(chromosome),
                                                                         gene_position = first(gene_position),
                                                                         N = first(N),
                                                                         phenotype_key = first(phenotype_key))
      
      trait2_sumstats = trait2_sumstats %>% group_by(gene) %>% summarise(betas = list(beta),
                                                                         variant_variances = list(variant_variance),
                                                                         chromosome = first(chromosome),
                                                                         gene_position = first(gene_position),
                                                                         N = first(N),
                                                                         phenotype_key = first(phenotype_key))
    }
    
  
  
  
  
  #compute w_t_beta, burden_score, and overdispersion from betas and variant_variances
  if (is.null(custom_weights)){
    if (null_stats) {
      trait1_sumstats$w_t_beta = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]*trait1_sumstats$betas[[x]] * sample(c(-1,1),length(trait1_sumstats$betas[[x]]),replace = TRUE)))
    } else {
      trait1_sumstats$w_t_beta = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]*trait1_sumstats$betas[[x]]))
    }
    trait1_sumstats$burden_score = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]))
    trait1_sumstats$overdispersion = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]^2)/sum(trait1_sumstats$variant_variances[[x]]))
    
    trait1_sumstats = trait1_sumstats[trait1_sumstats$burden_score > 0 & is.finite(trait1_sumstats$burden_score) & is.finite(trait1_sumstats$w_t_beta) & is.finite(trait1_sumstats$overdispersion),]
    
    if (!is.null(trait2_sumstats)){
      if (null_stats) {
        trait2_sumstats$w_t_beta = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]*trait2_sumstats$betas[[x]] * sample(c(-1,1),length(trait2_sumstats$betas[[x]]),replace = TRUE)))
      } else {
        trait2_sumstats$w_t_beta = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]*trait2_sumstats$betas[[x]]))
      }
      trait2_sumstats$burden_score = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]))
      trait2_sumstats$overdispersion = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]^2)/sum(trait2_sumstats$variant_variances[[x]]))
      trait2_sumstats = trait2_sumstats[trait2_sumstats$burden_score > 0 & is.finite(trait2_sumstats$burden_score) & is.finite(trait2_sumstats$w_t_beta) & is.finite(trait2_sumstats$overdispersion),]
      
    }
  } 
  
  
  else {
    trait1_sumstats$Wg_ug = sapply(1:nrow(trait1_sumstats), function(x) sqrt(trait1_sumstats$variant_variances[[x]]) * sqrt(trait1_sumstats$custom_weights[[x]]))
    
    if (null_stats) {
      trait1_sumstats$w_t_beta = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$Wg_ug[[x]]* sqrt(trait1_sumstats$variant_variances[[x]])*trait1_sumstats$betas[[x]]*sample(c(-1,1),length(trait1_sumstats$betas[[x]]), replace = TRUE)))
    } else {
      trait1_sumstats$w_t_beta = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$Wg_ug[[x]]* sqrt(trait1_sumstats$variant_variances[[x]])*trait1_sumstats$betas[[x]]))
    }
    trait1_sumstats$burden_score = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$Wg_ug[[x]]^2))
    trait1_sumstats$overdispersion = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]^2)/sum(trait1_sumstats$variant_variances[[x]]))
    
    trait1_sumstats = trait1_sumstats[trait1_sumstats$burden_score > 0 & is.finite(trait1_sumstats$burden_score) & is.finite(trait1_sumstats$w_t_beta) & is.finite(trait1_sumstats$overdispersion),]
    
    if (!is.null(trait2_sumstats)){
      trait2_sumstats$Wg_ug = sapply(1:nrow(trait2_sumstats), function(x) sqrt(trait2_sumstats$variant_variances[[x]]) * sqrt(trait2_sumstats$custom_weights[[x]]))
      if (null_stats) {
        trait2_sumstats$w_t_beta = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$Wg_ug[[x]]* sqrt(trait2_sumstats$variant_variances[[x]])*trait2_sumstats$betas[[x]]*sample(c(-1,1),length(trait2_sumstats$betas[[x]]), replace = TRUE)))
      } else {
        trait2_sumstats$w_t_beta = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$Wg_ug[[x]]* sqrt(trait2_sumstats$variant_variances[[x]])*trait2_sumstats$betas[[x]]))
      }
      trait2_sumstats$burden_score = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$Wg_ug[[x]]^2))
      trait2_sumstats$overdispersion = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]^2)/sum(trait2_sumstats$variant_variances[[x]]))

      trait2_sumstats = trait2_sumstats[trait2_sumstats$burden_score > 0 & is.finite(trait2_sumstats$burden_score) & is.finite(trait2_sumstats$w_t_beta) & is.finite(trait2_sumstats$overdispersion),]

    }
  }
  
  
  print("fitting model")
  
  
  
  if (mode == "univariate"){
    output = BHR_h2(trait1_sumstats, annotations, num_blocks, genomewide_correction, fixed_genes,output_jackknife_h2, overdispersion, all_models,num_null_conditions,slope_correction, gwc_exclusion, intercept)
    return(output)
  } else if (mode == "bivariate"){
    output = BHR_rg(trait1_sumstats = trait1_sumstats, trait2_sumstats = trait2_sumstats,annotations =  annotations, num_blocks = num_blocks,genomewide_correction = genomewide_correction,overdispersion = overdispersion,num_null_conditions = 0,output_jackknife_rg = output_jackknife_rg,fixed_genes = fixed_genes)
    return(output)
  } else if (mode == "aggregate"){
    output = BHR_meta(ss_list_trait1, trait_list, annotations, num_blocks, genomewide_correction, fixed_genes,  overdispersion, all_models, num_null_conditions, slope_correction, gwc_exclusion)
    return(output)
    } else if (mode == "aggregate-rg"){
      output = BHR_meta_rg(ss_list_trait1,ss_list_trait2, annotations, num_blocks)
      return(output)
    } else {
      return("Please enter a valid mode among: ['univariate','bivariate','aggregate', 'aggregate-rg']")
  }
}
