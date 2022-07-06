BHR <- function(mode = NULL,
                trait1_sumstats = NULL,
                trait2_sumstats = NULL,
                annotations, 
                num_blocks = 100, 
                genomewide_correction = FALSE, 
                fixed_genes = NULL,
                output_jackknife_h2 = FALSE,
                ss_list_trait1 = NULL,
                ss_list_trait2 = NULL,
                trait_list = NULL,
                overdispersion = FALSE,
                all_models = FALSE,
                slope_correction = NULL,
                num_null_conditions = 5,
                gwc_exclusion = NULL) {
  
  #compute w_t_beta, burden_score, and overdispersion from betas and variant_variances
  
  trait1_sumstats$w_t_beta = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]*trait1_sumstats$betas[[x]]))
  trait1_sumstats$burden_score = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]))
  trait1_sumstats$overdispersion = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]^2)/sum(trait1_sumstats$variant_variances[[x]]))
  
  if (!is.null(trait2_sumstats)){
    trait2_sumstats$w_t_beta = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]*trait2_sumstats$betas[[x]]))
    trait2_sumstats$burden_score = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]))
    trait2_sumstats$overdispersion = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]^2)/sum(trait2_sumstats$variant_variances[[x]]))
  }
  
  if (mode == "univariate"){
    output = BHR_h2(trait1_sumstats, annotations, num_blocks, genomewide_correction, fixed_genes,output_jackknife_h2, overdispersion, all_models,num_null_conditions,slope_correction, gwc_exclusion)
    return(output)
  } else if (mode == "bivariate"){
    output = BHR_rg(trait1_sumstats = trait1_sumstats, trait2_sumstats = trait2_sumstats,annotations =  annotations, num_blocks = num_blocks,genomewide_correction = genomewide_correction,overdispersion = overdispersion,num_null_conditions = 0,output_jackknife_rg = FALSE,fixed_genes = fixed_genes)
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
