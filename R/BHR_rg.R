
BHR_rg <- function(trait1_sumstats, 
                   trait2_sumstats, 
                   annotations, 
                   num_blocks = 100,
                   genomewide_correction = TRUE,
                   overdispersion, 
                   num_null_conditions = 0,
                   use_null_conditions_rg = FALSE,
                   output_jackknife_rg = FALSE, 
                   fixed_genes = NULL,
                   intercept = intercept,
                   rg_random_se_estimator = "jackknife",
                   log,
                   start_time) {
  trait1_sumstats = trait1_sumstats[!is.na(trait1_sumstats$chromosome),]
  trait2_sumstats = trait2_sumstats[!is.na(trait2_sumstats$chromosome),]

  consensus_genes = intersect(trait1_sumstats$gene, trait2_sumstats$gene)
  trait1_sumstats = trait1_sumstats[match(consensus_genes,trait1_sumstats$gene),]
  trait2_sumstats = trait2_sumstats[match(consensus_genes,trait2_sumstats$gene),]

  trait1_sumstats$gamma_sq = trait1_sumstats$w_t_beta^2/trait1_sumstats$burden_score
  trait2_sumstats$gamma_sq = trait2_sumstats$w_t_beta^2/trait2_sumstats$burden_score

  trait1_sumstats$gamma = trait1_sumstats$w_t_beta/sqrt(trait1_sumstats$burden_score)
  trait2_sumstats$gamma = trait2_sumstats$w_t_beta/sqrt(trait2_sumstats$burden_score)

  trait1_sumstats$burden_score_sqrt = sqrt(trait1_sumstats$burden_score)
  trait2_sumstats$burden_score_sqrt = sqrt(trait2_sumstats$burden_score)


  merged_annotations <- annotations %>% reduce(inner_join, by = "gene")

  if (is.null(fixed_genes)){
    trait1_sig_genes = trait1_sumstats$gene[(trait1_sumstats$N - 1) * (trait1_sumstats$gamma_sq) > qchisq(p = 1 - (0.05/nrow(trait1_sumstats)), df = 1)]
    trait2_sig_genes = trait2_sumstats$gene[(trait2_sumstats$N - 1) * (trait2_sumstats$gamma_sq) > qchisq(p = 1 - (0.05/nrow(trait2_sumstats)), df = 1)]
    sig_genes = union(trait1_sig_genes,trait2_sig_genes)
  } else{
    sig_genes = fixed_genes
  }

  if (use_null_conditions_rg) {
    heritability_trait1 <- BHR_h2(trait1_sumstats, annotations = annotations, num_blocks = num_blocks, fixed_genes = sig_genes, genomewide_correction = FALSE, output_jackknife_h2 = TRUE, overdispersion = overdispersion, num_null_conditions = num_null_conditions, slope_correction = FALSE, all_models = TRUE, gwc_exclusion = NULL, intercept = intercept, start_time = Sys.time(), log = log)
    heritability_trait2 <- BHR_h2(trait2_sumstats, annotations = annotations, num_blocks = num_blocks, fixed_genes = sig_genes, genomewide_correction = FALSE,output_jackknife_h2 = TRUE, overdispersion = overdispersion, num_null_conditions = num_null_conditions, slope_correction = FALSE, all_models = TRUE, gwc_exclusion = NULL, intercept = intercept, start_time = Sys.time(), log = log)
  } else { 
    heritability_trait1 <- BHR_h2(trait1_sumstats, annotations = annotations, num_blocks = num_blocks, fixed_genes = sig_genes, genomewide_correction = FALSE, output_jackknife_h2 = TRUE, overdispersion = overdispersion, num_null_conditions = 0, slope_correction = FALSE, all_models = TRUE, gwc_exclusion = NULL, intercept = intercept, start_time = Sys.time(), log = log)
    heritability_trait2 <- BHR_h2(trait2_sumstats, annotations = annotations, num_blocks = num_blocks, fixed_genes = sig_genes, genomewide_correction = FALSE,output_jackknife_h2 = TRUE, overdispersion = overdispersion, num_null_conditions = 0, slope_correction = FALSE, all_models = TRUE, gwc_exclusion = NULL, intercept = intercept, start_time = Sys.time(), log = log)
    }
  #add null moment conditions
  trait1_sumstats = trait1_sumstats[,c("gene","N","gamma_sq","gamma","w_t_beta","burden_score","burden_score_sqrt","overdispersion","chromosome","gene_position")]
  trait1_sumstats$true = TRUE

  trait2_sumstats = trait2_sumstats[,c("gene","N","gamma_sq","gamma","w_t_beta","burden_score","burden_score_sqrt","overdispersion","chromosome","gene_position")]
  trait2_sumstats$true = TRUE

  filter = trait1_sumstats$burden_score > 0 & trait2_sumstats$burden_score > 0
  trait1_sumstats = trait1_sumstats[filter,]
  trait2_sumstats = trait2_sumstats[filter,]


  pair_sumstats = data.frame(gene = trait1_sumstats$gene,
                             w_t_beta = trait1_sumstats$w_t_beta * trait2_sumstats$w_t_beta,
                             gamma_sq = trait1_sumstats$gamma*trait2_sumstats$gamma,
                             burden_score = trait1_sumstats$burden_score_sqrt * trait2_sumstats$burden_score_sqrt,
                             overdispersion = sqrt(trait1_sumstats$overdispersion) * sqrt(trait2_sumstats$overdispersion),
                             chromosome = trait1_sumstats$chromosome,
                             gene_position = trait1_sumstats$gene_position,
                             true = trait1_sumstats$true)


  pair_sumstats <- merge(pair_sumstats, merged_annotations, by.x = "gene", by.y = "gene")
  pair_sumstats = pair_sumstats[with(pair_sumstats, order(chromosome, gene_position)),]
  block = ceiling((1:nrow(pair_sumstats[!(pair_sumstats$gene %in% sig_genes),]))/(nrow(pair_sumstats[!(pair_sumstats$gene %in% sig_genes),])/num_blocks))
  genetic_covariance_random <- randomeffects_jackknife(pair_sumstats[!(pair_sumstats$gene %in% sig_genes),],merged_annotations,num_blocks,block,genomewide_correction = genomewide_correction, output_jackknife_h2 = TRUE, overdispersion = overdispersion,slope_correction = FALSE, bivariate = TRUE, intercept = intercept)

  pair_sumstats_true = pair_sumstats[pair_sumstats$true,]
  trait1_sumstats_true = trait1_sumstats[trait1_sumstats$true,]
  trait2_sumstats_true = trait2_sumstats[trait2_sumstats$true,]

  block_true = ceiling((1:nrow(pair_sumstats_true))/(nrow(pair_sumstats_true)/num_blocks))


  #random effects model results
  #point estimate
  rho_random = genetic_covariance_random$heritabilities[1,ncol(genetic_covariance_random$heritabilities)]
  rho_random_se = genetic_covariance_random$heritabilities[2,ncol(genetic_covariance_random$heritabilities)]

  h2_trait1_random = heritability_trait1$subthreshold_genes$heritabilities[1,ncol(heritability_trait1$subthreshold_genes$heritabilities)]
  h2_trait2_random = heritability_trait2$subthreshold_genes$heritabilities[1,ncol(heritability_trait2$subthreshold_genes$heritabilities)]
  rg_random = rho_random/sqrt(h2_trait1_random * h2_trait2_random)

  #se
  rho_jackknife = genetic_covariance_random$jackknife_h2[nrow(genetic_covariance_random$jackknife_h2),]
  h2_trait1_jackknife = heritability_trait1$subthreshold_genes$jackknife_h2[nrow(heritability_trait1$subthreshold_genes$jackknife_h2),]
  h2_trait2_jackknife = heritability_trait2$subthreshold_genes$jackknife_h2[nrow(heritability_trait2$subthreshold_genes$jackknife_h2),]
  
  if (rg_random_se_estimator == "delta") {
  gradient_rg_subthreshold = c(1/sqrt(h2_trait1_random*h2_trait2_random),
                  (-rho_random*h2_trait2_random)/(2 * (h2_trait1_random*h2_trait2_random)^(-3/2)),
                  (-rho_random*h2_trait1_random)/(2 * (h2_trait1_random*h2_trait2_random)^(-3/2)))
  
  jackknife_rg_subthreshold_mat = cbind(rho_jackknife - mean(rho_jackknife), h2_trait1_jackknife - mean(h2_trait1_jackknife),h2_trait2_jackknife - mean(h2_trait2_jackknife))
  S_rg_subthreshold = ((num_blocks-1)/num_blocks)*(t(jackknife_rg_subthreshold_mat) %*% jackknife_rg_subthreshold_mat)
  
  rg_random_var = t(gradient_rg_subthreshold) %*% S_rg_subthreshold %*% gradient_rg_subthreshold
  rg_random_se = sqrt(rg_random_var)
  
  } else if (rg_random_se_estimator == "jackknife") {
    rg_jackknife = rho_jackknife/sqrt(h2_trait1_jackknife * h2_trait2_jackknife)
    rg_random_se = sqrt(((num_blocks-1)/num_blocks) * sum((rg_jackknife - mean(rg_jackknife))^2))
    
  }
  


  if (length(sig_genes) > 0) {
    #mixed effects model
    #point estimate
    pair_sumstats_sig = pair_sumstats_true[pair_sumstats_true$gene %in% sig_genes,]
    trait1_sumstats_sig = trait1_sumstats_true[trait1_sumstats_true$gene %in% sig_genes,]
    trait2_sumstats_sig = trait2_sumstats_true[trait2_sumstats_true$gene %in% sig_genes,]

    #FIX: for each gene, divide by product of norms of w
    rho_mixed = rho_random + sum(pair_sumstats_sig$gamma_sq)
    trait1_h2_mixed = h2_trait1_random + sum(trait1_sumstats_sig$w_t_beta^2/trait1_sumstats_sig$burden_score)
    trait2_h2_mixed = h2_trait2_random + sum(trait2_sumstats_sig$w_t_beta^2/trait2_sumstats_sig$burden_score)

    rg_mixed = rho_mixed/sqrt(trait1_h2_mixed*trait2_h2_mixed)

    #standard error
    #gradient of the genetic correlation
    #first, dg_duk terms
    u = trait1_sumstats_sig$w_t_beta/sqrt(trait1_sumstats_sig$burden_score)
    v = trait2_sumstats_sig$w_t_beta/sqrt(trait2_sumstats_sig$burden_score)
    h11 = h2_trait1_random
    h22 = h2_trait2_random
    h12 = rho_random
    dg_duks <- sapply(1:length(sig_genes), function(x) {
      u_k = u[x]
      v_k = v[x]

      partial = ((((t(u) %*% v)+h12) * (-u_k))/(sqrt((t(v) %*% v)+h22)*(((t(u) %*% u)+h11)^(3/2)))) + ((v_k)/((sqrt((t(v) %*% v)+h22))*(sqrt((t(u) %*% u)+h11))))

      return(partial)
    })

    dg_dvks <- sapply(1:length(sig_genes), function(x) {
      u_k = u[x]
      v_k = v[x]

      partial = ((((t(v) %*% u)+h12) * (-v_k))/(sqrt((t(u) %*% u)+h11)*(((t(v) %*% v)+h22)^(3/2)))) + ((u_k)/((sqrt((t(u) %*% u)+h11))*(sqrt((t(v) %*% v)+h22))))

      return(partial)
    })

    dg_dh12 = 1/((sqrt((t(u) %*% u)+h11))*(sqrt((t(v) %*% v)+h22)))
    dg_dh11 = (((t(v) %*% u)+h12)/sqrt((t(v) %*% v)+h22))*(-0.5)*(1/(((t(u) %*% u)+h11)^(3/2)))
    dg_dh22 = (((t(u) %*% v)+h12)/sqrt((t(u) %*% u)+h11))*(-0.5)*(1/(((t(v) %*% v)+h22)^(3/2)))

    gradient = c(dg_duks,dg_dvks,dg_dh12,dg_dh11,dg_dh22)

    #now, calculate covariance matrix of arguments of genetic correlation

    trait1_covariance = heritability_trait1$subthreshold_genes$intercept * diag(length(sig_genes))
    trait2_covariance = heritability_trait2$subthreshold_genes$intercept * diag(length(sig_genes))
    crosstrait_covariance = genetic_covariance_random$intercept * diag(length(sig_genes))

    fixed_covariance = rbind(cbind(trait1_covariance,crosstrait_covariance), cbind(crosstrait_covariance,trait2_covariance))

    jackknife_heritabilities_covariance = cbind(rho_jackknife,h2_trait1_jackknife,h2_trait2_jackknife)
    jackknife_heritabilities_covariance <- sapply(1:ncol(jackknife_heritabilities_covariance), function(x) jackknife_heritabilities_covariance[,x] - mean(jackknife_heritabilities_covariance[,x]))

    jackknife_heritabilities_covariancematrix = ((num_blocks-1)/num_blocks)*(t(jackknife_heritabilities_covariance) %*% jackknife_heritabilities_covariance)

    S <- matrix(data = 0, nrow = nrow(fixed_covariance)+nrow(jackknife_heritabilities_covariancematrix), ncol = nrow(fixed_covariance)+nrow(jackknife_heritabilities_covariancematrix))
    S[1:nrow(fixed_covariance),1:nrow(fixed_covariance)] <- fixed_covariance
    S[(nrow(S) - nrow(jackknife_heritabilities_covariancematrix)+1):nrow(S),(nrow(S) - nrow(jackknife_heritabilities_covariancematrix)+1):nrow(S)] <- jackknife_heritabilities_covariancematrix


    rg_mixed_se = sqrt(t(gradient) %*% S %*% gradient)

    rho_mixed_se = sqrt(jackknife_heritabilities_covariancematrix[1,1] + (genetic_covariance_random$intercept * length(sig_genes)))

  } else {
    rg_mixed = rg_random
    rg_mixed_se = rg_random_se
    rho_mixed = rho_random
    rho_mixed_se = rho_random_se
    #   fixed_results = NA
  }
  if (output_jackknife_rg){
    output = list(trait1 = heritability_trait1,
                trait2 = heritability_trait2,
                rg = list(rg_subthreshold = rg_random,
                          rg_subthreshold_se = rg_random_se,
                          rg_mixed = rg_mixed,
                          rg_mixed_se = rg_mixed_se,
                          rho_mixed = rho_mixed,
                          rho_mixed_se = rho_mixed_se,
                          rho_subthreshold = rho_random,
                          rho_subthreshold_se = rho_random_se,
                          rho_jackknife = rho_jackknife,
                          rg_jackknife = rg_jackknife,
                          intercept = genetic_covariance_random$intercept,
                          sig_genes = if (length(sig_genes) > 0) {pair_sumstats_sig$gene} else {NA},
                          #         fractions_mixed = fixed_results$fractions,
                          fractions_subthreshold = genetic_covariance_random$fractions,
                          enrichments_subthreshold = genetic_covariance_random$enrichments_final),
                run_time = Sys.time() - start_time)
  } else {
    output = list(trait1 = heritability_trait1,
                  trait2 = heritability_trait2,
                  rg = list(rg_subthreshold = rg_random,
                            rg_subthreshold_se = rg_random_se,
                            rg_mixed = rg_mixed,
                            rg_mixed_se = rg_mixed_se,
                            rho_subthreshold = rho_random,
                            rho_subthreshold_se = rho_random_se,
                            rho_mixed = rho_mixed,
                            rho_mixed_se = rho_mixed_se,
                            intercept = genetic_covariance_random$intercept,
                            sig_genes = if (length(sig_genes) > 0) {pair_sumstats_sig$gene} else {NA},
                            #         fractions_mixed = fixed_results$fractions,
                            fractions_subthreshold = genetic_covariance_random$fractions,
                            enrichments_subthreshold = genetic_covariance_random$enrichments_final),
                  run_time = Sys.time() - start_time)
  }


  return(output)
}
