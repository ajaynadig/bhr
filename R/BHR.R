BHR <- function(mode = NULL,
                trait1_sumstats = NULL,
                trait2_sumstats = NULL,
                annotations,
                num_blocks = 100,
                genomewide_correction = FALSE,
                gwc_exclusion = NULL,
                fixed_genes = NULL,
                output_jackknife_h2 = FALSE,
                output_jackknife_rg = FALSE,
                ss_list = NULL,
                trait_list = NULL,
                overdispersion = FALSE,
                all_models = FALSE,
                slope_correction = NULL,
                num_null_conditions = 5,
                use_null_conditions_rg = FALSE,
                custom_weights = FALSE,
                null_stats = FALSE,
                intercept = TRUE,
                custom_variant_variances = FALSE,
                rg_random_se_estimator = "jackknife") {
  start_time = Sys.time()
  message("Burden Heritability Regression\nDaniel Weiner, Ajay Nadig, and Luke O'Connor, 2023")
  message(paste0("Running BHR at ",start_time))
  
  if (mode == "aggregate"){
    output = BHR_meta(ss_list, 
                      trait_list, 
                      annotations, 
                      num_blocks, 
                      genomewide_correction, 
                      fixed_genes,  
                      overdispersion, 
                      all_models, 
                      num_null_conditions, 
                      slope_correction, 
                      gwc_exclusion,
                      intercept,
                      custom_variant_variances = custom_variant_variances,
                      start_time = start_time)
    message(paste0("BHR finished at ",Sys.time()))
    return(output)
  } else {
  
  #Calculate variance variances from AF, if not provided
  if (custom_variant_variances == FALSE){
    if (is.null(trait2_sumstats)){
      trait1_sumstats$variant_variance = 2*trait1_sumstats$AF*(1-trait1_sumstats$AF)
    } else {
      trait1_sumstats$variant_variance = 2*trait1_sumstats$AF*(1-trait1_sumstats$AF)
      trait2_sumstats$variant_variance = 2*trait2_sumstats$AF*(1-trait2_sumstats$AF)
    }
  } 
  
  #Check range of allele frequencies
  allele_frequencies_trait1 = (2  - sqrt(4 - 8*trait1_sumstats$variant_variance))/4
  allele_frequencies_trait2 = (2  - sqrt(4 - 8*trait2_sumstats$variant_variance))/4
  
  max_frequency_trait1 = max(allele_frequencies_trait1)
  min_frequency_trait1 = min(allele_frequencies_trait1)
  diff_log10_trait1 = log10(max_frequency_trait1) - log10(min_frequency_trait1)
  
  message(paste("For trait 1, MAF ranges from", as.character(signif(min_frequency_trait1,2)),"to", as.character(signif(max_frequency_trait1,2))))
  
  if (diff_log10_trait1 <= 2){
    message("...seems reasonable")
  } else {
    message("MAF ranges over more than two orders of magnitude for trait 1. We recommend splitting into finer MAF bins. See documentation")
  }
  
  if (!is.null(trait2_sumstats)){
    max_frequency_trait2 = max(allele_frequencies_trait2)
    min_frequency_trait2 = min(allele_frequencies_trait2)
    diff_log10_trait2 = log10(max_frequency_trait2) - log10(min_frequency_trait2)
    
    message(paste("For trait 2, MAF ranges from", as.character(signif(min_frequency_trait2,2)),"to", as.character(signif(max_frequency_trait2,2))))
    
    if (diff_log10_trait2 <=2){
      message("...seems reasonable")
    } else {
      message("MAF ranges over more than two of magnitude for trait 2. We recommend splitting into finer MAF bins. See documentation")
    }
  }
  
  
  #Aggregate variant-row summary statistics into gene-row summary statistics
    if (is.null(trait2_sumstats)){
      trait1_sumstats = trait1_sumstats %>% group_by(gene) %>% summarise(betas = list(beta),
                                                                         variant_variances = list(variant_variance),
                                                                         chromosome = first(chromosome),
                                                                         gene_position = first(gene_position),
                                                                         N = first(N),
                                                                         phenotype_key = first(phenotype_key),
                                                                         num_variants = length(variant_variance))
    } else{
      trait1_sumstats = trait1_sumstats %>% group_by(gene) %>% summarise(betas = list(beta),
                                                                         variant_variances = list(variant_variance),
                                                                         chromosome = first(chromosome),
                                                                         gene_position = first(gene_position),
                                                                         N = first(N),
                                                                         phenotype_key = first(phenotype_key),
                                                                         num_variants = length(variant_variance))
      
      trait2_sumstats = trait2_sumstats %>% group_by(gene) %>% summarise(betas = list(beta),
                                                                         variant_variances = list(variant_variance),
                                                                         chromosome = first(chromosome),
                                                                         gene_position = first(gene_position),
                                                                         N = first(N),
                                                                         phenotype_key = first(phenotype_key),
                                                                         num_variants = length(variant_variance))
    }
    
  
  #compute w_t_beta, burden_score, and overdispersion from betas and variant_variances
  if (!custom_weights){
    if (null_stats) {
      trait1_sumstats$w_t_beta = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]*trait1_sumstats$betas[[x]] * sample(c(-1,1),length(trait1_sumstats$betas[[x]]),replace = TRUE)))
    } else {
      trait1_sumstats$w_t_beta = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]*trait1_sumstats$betas[[x]]))
    }
    trait1_sumstats$burden_score = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]))
    trait1_sumstats$cumulative_variance = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]))
    trait1_sumstats$overdispersion = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]^2)/sum(trait1_sumstats$variant_variances[[x]]))
    
    trait1_sumstats = trait1_sumstats[trait1_sumstats$burden_score > 0 & is.finite(trait1_sumstats$burden_score) & is.finite(trait1_sumstats$w_t_beta) & is.finite(trait1_sumstats$overdispersion),]
    
    if (!is.null(trait2_sumstats)){
      if (null_stats) {
        trait2_sumstats$w_t_beta = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]*trait2_sumstats$betas[[x]] * sample(c(-1,1),length(trait2_sumstats$betas[[x]]),replace = TRUE)))
      } else {
        trait2_sumstats$w_t_beta = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]*trait2_sumstats$betas[[x]]))
      }
      trait2_sumstats$burden_score = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]))
      trait2_sumstats$cumulative_variance = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]))
      trait2_sumstats$overdispersion = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]^2)/sum(trait2_sumstats$variant_variances[[x]]))
      trait2_sumstats = trait2_sumstats[trait2_sumstats$burden_score > 0 & is.finite(trait2_sumstats$burden_score) & is.finite(trait2_sumstats$w_t_beta) & is.finite(trait2_sumstats$overdispersion),]
      
    }
  }  else {
    trait1_sumstats$Wg_ug = sapply(1:nrow(trait1_sumstats), function(x) sqrt(trait1_sumstats$variant_variances[[x]]) * sqrt(trait1_sumstats$custom_weights[[x]]))
    
    if (null_stats) {
      trait1_sumstats$w_t_beta = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$Wg_ug[[x]]* sqrt(trait1_sumstats$variant_variances[[x]])*trait1_sumstats$betas[[x]]*sample(c(-1,1),length(trait1_sumstats$betas[[x]]), replace = TRUE)))
    } else {
      trait1_sumstats$w_t_beta = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$Wg_ug[[x]]* sqrt(trait1_sumstats$variant_variances[[x]])*trait1_sumstats$betas[[x]]))
    }
    trait1_sumstats$burden_score = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$Wg_ug[[x]]^2))
    trait1_sumstats$cumulative_variance = sapply(1:nrow(trait1_sumstats), function(x) sum(trait1_sumstats$variant_variances[[x]]))
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
      trait2_sumstats$cumulative_variance = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]))
      trait2_sumstats$overdispersion = sapply(1:nrow(trait2_sumstats), function(x) sum(trait2_sumstats$variant_variances[[x]]^2)/sum(trait2_sumstats$variant_variances[[x]]))

      trait2_sumstats = trait2_sumstats[trait2_sumstats$burden_score > 0 & is.finite(trait2_sumstats$burden_score) & is.finite(trait2_sumstats$w_t_beta) & is.finite(trait2_sumstats$overdispersion),]

    }
  }

  #create log
  log = list(mode = mode,
             num_blocks = num_blocks,
             genomewide_correction = genomewide_correction,
             gwc_exclusion = gwc_exclusion,
             fixed_genes = fixed_genes,
             output_jackknife_h2 = output_jackknife_h2,
             output_jackknife_rg = output_jackknife_rg,
             overdispersion = overdispersion,
             all_models = all_models,
             slope_correction = slope_correction,
             num_null_conditions = num_null_conditions,
             custom_weights = custom_weights,
             null_stats = null_stats,
             intercept = intercept,
             custom_variant_variances = custom_variant_variances)
  
  if (mode == "univariate"){
    output = BHR_h2(trait1_sumstats, 
                    annotations, 
                    num_blocks, 
                    genomewide_correction, 
                    fixed_genes,
                    output_jackknife_h2, 
                    overdispersion, 
                    all_models,
                    num_null_conditions,
                    slope_correction, 
                    gwc_exclusion, 
                    intercept,
                    log,
                    start_time)
    message(paste0("BHR finished at ",Sys.time()))
    return(output)
  } else if (mode == "bivariate"){
    output = BHR_rg(trait1_sumstats = trait1_sumstats, 
                    trait2_sumstats = trait2_sumstats,
                    annotations =  annotations, 
                    num_blocks = num_blocks,
                    genomewide_correction = genomewide_correction,
                    overdispersion = overdispersion,
                    num_null_conditions = num_null_conditions,
                    use_null_conditions_rg = use_null_conditions_rg,
                    output_jackknife_rg = output_jackknife_rg,
                    fixed_genes = fixed_genes, 
                    rg_random_se_estimator = rg_random_se_estimator,
                    log = log,
                    intercept = intercept,
                    start_time = start_time)
    message(paste0("BHR finished at ",Sys.time()))
    return(output)
    } else {
      return("Please enter a valid mode among: ['univariate','bivariate','aggregate']")
      message(paste0("BHR finished at ",Sys.time()))
  }
  }
}



