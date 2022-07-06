library(tidyverse)

BHR_h2 <- function(sumstats, annotations, num_blocks, genomewide_correction, fixed_genes, output_jackknife_h2, overdispersion, all_models, num_null_conditions,slope_correction, gwc_exclusion) {
  set.seed(363)
  #Step 1: Wrangling, functions, and other code relevant to both univariate and rg cases.

  if (genomewide_correction == TRUE & is.null(gwc_exclusion)){
    print("Running GWC based on all genes")

    mu_genome = sum(sumstats$w_t_beta)/sum(sumstats$burden_score)
    per_gene_correction = mu_genome * sqrt(sumstats$burden_score)
    sumstats$gamma_sq = ((sumstats$w_t_beta / sqrt(sumstats$burden_score)) - per_gene_correction)^2

    print(paste0("Mean genome-wide per allele effect size: ", signif(mu_genome, 2)))

  } else if (genomewide_correction == TRUE & !is.null(gwc_exclusion)) {
    print(paste0("Running GWC based on all genes, except in annotation ", gwc_exclusion))

    merged_annotations <- annotations %>% reduce(inner_join, by = "gene")
    exclusion_genes = merged_annotations[merged_annotations[,gwc_exclusion] == 1,"gene"]
    print(paste0(length(exclusion_genes), " genes excluded from annotation"))
    sumstats_edited = sumstats[!(sumstats$gene %in% exclusion_genes),]

    mu_genome = sum(sumstats_edited$w_t_beta)/sum(sumstats_edited$burden_score)
    per_gene_correction = mu_genome * sqrt(sumstats$burden_score)
    sumstats$gamma_sq = ((sumstats$w_t_beta / sqrt(sumstats$burden_score)) - per_gene_correction)^2

    print(paste0("Mean genome-wide per allele effect size: ", signif(mu_genome, 2)))

  } else  {
    mu_genome = sum(sumstats$w_t_beta)/sum(sumstats$burden_score)
    sumstats$gamma_sq = sumstats$w_t_beta^2/sumstats$burden_score
  }

  if (num_null_conditions >0) {
    nullstats = sumstats[,c("gene","N","overdispersion","chromosome","gene_position","variant_variances","betas","trait_type")]
    truestats = sumstats[,c("gene","N","gamma_sq","w_t_beta","burden_score","overdispersion","chromosome","gene_position","trait_type")]
    truestats$true = TRUE
  } else{
    truestats = sumstats[,c("gene","N","gamma_sq","w_t_beta","burden_score","overdispersion","chromosome","gene_position")]
    truestats$true = TRUE
  }

  sumstats = truestats
  if (num_null_conditions > 0){
    for (nullcondition in 1:num_null_conditions){
      w_true = sapply(1:length(nullstats$variant_variances), function(x) sqrt(nullstats$variant_variances[[x]]))

      w_null <- sapply(1:length(nullstats$variant_variances), function(x){
        w_true[[x]]*sample(c(-1,1),length(w_true[[x]]), replace = TRUE)})


      get_deltasq <- function(row){
        beta_tilde = nullstats$betas[[row]]
        beta = beta_tilde * w_true[[row]]
        # if (nullstats$trait_type[row] == "continuous"){
        #   beta_tilde = nullstats$betas[[row]]
        #   beta = beta_tilde * w_true[[row]]
        # } else if (nullstats$trait_type[row] == "categorical"){
        #   beta = nullstats$betas[[row]]
        # }
        return((sum(w_null[[row]]*beta)^2)/sum(w_null[[row]]*w_null[[row]]))
      }

      nullstats$gamma_sq <- sapply(1:nrow(nullstats),get_deltasq)

      get_null_w_t_beta <- function(row){
        beta_tilde = nullstats$betas[[row]]
        beta = beta_tilde * w_true[[row]]
        # if (nullstats$trait_type[row] == "continuous"){
        #   beta_tilde = nullstats$betas[[row]]
        #   beta = beta_tilde * w_true[[row]]
        # } else if (nullstats$trait_type[row] == "categorical"){
        #   beta = nullstats$betas[[row]]
        # }
        return((sum(w_null[[row]]*beta)))
      }

      nullstats$w_t_beta <- sapply(1:nrow(nullstats),get_null_w_t_beta)

      get_null_burdenscore<- function(row){
        return((sum(w_null[[row]]*w_true[[row]])^2)/sum(w_null[[row]]*w_null[[row]]))
      }

      nullstats$burden_score <- sapply(1:nrow(nullstats),get_null_burdenscore)
      nullstats_out = nullstats[,c("gene","N","gamma_sq","w_t_beta","burden_score","overdispersion","chromosome","gene_position","trait_type")]
      nullstats_out$true = FALSE

      sumstats = rbind(sumstats,nullstats_out)
    }
  }

  #wrangle annotations
  merged_annotations <- annotations %>% reduce(inner_join, by = "gene")

  #Step 2: Calculate heritability of trait 1, and if univariate, return results.
  sumstats <- merge(sumstats, merged_annotations, by.x = "gene", by.y = "gene")
  sumstats = sumstats[with(sumstats, order(chromosome, gene_position)),]

  if (!is.null(fixed_genes)) {
    sig_chisq = sumstats$gene[sumstats$true] %in% fixed_genes
  } else {
    sig_chisq = (sumstats$N[sumstats$true] - 1) * (sumstats$gamma_sq[sumstats$true]) > qchisq(p = 1 - (0.05/nrow(sumstats[sumstats$true,])), df = 1)
  }

  sumstats_sig = sumstats[sumstats$true,][sig_chisq,]

  #Estimate h2 and h2_SE using only genes below GW threshold
  sumstats_sub = sumstats[!(sumstats$gene %in% sumstats_sig$gene),]
  block = ceiling((1:nrow(sumstats_sub))/(nrow(sumstats_sub)/num_blocks))

  subthreshold_genes_h2 <- randomeffects_jackknife(sumstats_sub, merged_annotations,num_blocks,block,genomewide_correction,output_jackknife_h2 = TRUE, overdispersion, slope_correction)

  #now that the regression is complete, we can get rid of the null burden moment conditions

  sumstats_true = sumstats[sumstats$true,]
  block_true = ceiling((1:nrow(sumstats_true))/(nrow(sumstats_true)/num_blocks))

  #Estimate h2 and h2_SE using the mixed model
  if (nrow(sumstats_sig) > 0){

    fixed_annotation_h2 = (sumstats_sig$gamma_sq) %*% as.matrix(sumstats_sig[,(ncol(sumstats_sig)-ncol(merged_annotations)+2):ncol(sumstats_sig)])
    fixed_h2_out <- c(fixed_annotation_h2, sum(sumstats_sig$gamma_sq))

    fixed_annotation_se = ((4*(sumstats_sig$gamma_sq))/sumstats_sig$N) %*% as.matrix(sumstats_sig[,(ncol(sumstats_sig)-ncol(merged_annotations)+2):ncol(sumstats_sig)])
    fixed_var_out <- c(fixed_annotation_se, sum((4*(sumstats_sig$gamma_sq))/sumstats_sig$N)) + (4/sumstats_sig$N[1]^2)

    sampling_noise_correction_factor = c(as.numeric((colSums((as.matrix(sumstats_sig[,(ncol(sumstats_sig)-ncol(merged_annotations)+2):ncol(sumstats_sig)]))))), nrow(sumstats_sig))/sumstats_true$N[1]
    mixed_h2 = subthreshold_genes_h2$heritabilities[1,1:ncol(merged_annotations)] + fixed_h2_out - sampling_noise_correction_factor
    mixed_h2_se = sqrt(subthreshold_genes_h2$heritabilities[2,1:ncol(merged_annotations)]^2 + fixed_var_out)

    h2_output_table <- rbind(mixed_h2,mixed_h2_se)
    colnames(h2_output_table)[ncol(h2_output_table)] = "total"
    rownames(h2_output_table) = c("h2","h2_se")

    #calculate mixed model fractions of heritability

    #first, get the point estimates

    mixed_model_fractionh2s = mixed_h2[1:(ncol(merged_annotations)-1)]/mixed_h2[ncol(merged_annotations)]

    #then, use delta method to get SEs
    #first, create covariance matrix of fixed effects and random effects
    #covariance matrix of fixed effects
    fixed_covariance = diag(nrow(sumstats_sig)) * subthreshold_genes_h2$intercept

    get_mixed_fraction_SEs <- function(annot_index){
      jackknife_fractions <- t(subthreshold_genes_h2$jackknife_h2[c(annot_index,nrow(subthreshold_genes_h2$jackknife_h2)),])
      jackknife_fractions <- sapply(1:ncol(jackknife_fractions), function(x) jackknife_fractions[,x] - mean(jackknife_fractions[,x]))
      mixed_fraction_covariance = ((num_blocks-1)/num_blocks)*(t(jackknife_fractions) %*% jackknife_fractions)

      S <- matrix(data = 0, nrow = nrow(fixed_covariance)+nrow(mixed_fraction_covariance), ncol = nrow(fixed_covariance)+nrow(mixed_fraction_covariance))
      S[1:nrow(fixed_covariance),1:nrow(fixed_covariance)] <- fixed_covariance
      S[(nrow(S) - nrow(mixed_fraction_covariance)+1):nrow(S),(nrow(S) - nrow(mixed_fraction_covariance)+1):nrow(S)] <- mixed_fraction_covariance

      #compute the gradient
      u = sumstats_sig$gamma_sq
      a = sumstats_sig[,(ncol(sumstats_true) - ncol(merged_annotations) + annot_index)]
      h_all = subthreshold_genes_h2$heritabilities[1,ncol(merged_annotations)]
      h_annot = subthreshold_genes_h2$heritabilities[1,annot_index]

      A = diag(length(a))
      diag(A) = a

      dg_duk_terms <- sapply(1:length(u), function(k) ((2*u[k]*a[k])/((t(u)%*%u) + h_all)) - ((2*u[k]*((t(u)%*%A%*%u) + h_annot))/((t(u)%*%u) + h_all)^2))

      dg_dhannot <- 1/((t(u)%*%u) + h_all)
      dg_dhall <- -((t(u)%*%A%*%u) + h_annot)/(((t(u)%*%u) + h_all)^2)

      gradient_fraction_annot = c(dg_duk_terms,dg_dhannot,dg_dhall)

      fraction_se <- t(gradient_fraction_annot) %*% S %*% gradient_fraction_annot

    }
    mixed_fraction_SEs = sqrt(sapply(1:(ncol(merged_annotations)-1), get_mixed_fraction_SEs))

    fixed_fraction_output_table <- rbind(mixed_model_fractionh2s,mixed_fraction_SEs)
    rownames(fixed_fraction_output_table) <- c("fraction_h2","fraction_h2_se")

    fraction_burden_score <- (sumstats_true$burden_score %*% data.matrix(sumstats_true[colnames(merged_annotations[-1])]))/sum(sumstats_true$burden_score)
    enrichments <- data.frame(sweep(fixed_fraction_output_table, 2, fraction_burden_score, `/`))
    rownames(enrichments) <- c('enrichments', 'enrichment_se')

    mixed_results <- list(heritabilities = h2_output_table,
                          fractions = fixed_fraction_output_table,
                          fraction_burden_score = fraction_burden_score,
                          enrichments = enrichments)
    number_genomesig_genes <- nrow(sumstats_sig)
    genomesig_genes <- as.character(sumstats_sig$gene)

    #get proportions of heritabillity explained by significant genes
    #point estimate

    #standard error through delta method
    u = sqrt((sumstats_sig$gamma_sq))
    u_t_u = t(u) %*% u
    h = subthreshold_genes_h2$heritabilities[1,ncol(merged_annotations)]
    #FIX
    frac_sig = u_t_u/mixed_h2[ncol(merged_annotations)]

    dg_duks <- sapply(1:nrow(sumstats_sig), function(k) {
      partial = ((-u_t_u * 2*u[k])/(u_t_u + h)^2) + ((2*u[k])/(u_t_u + h))
      return(partial)
    })
    dg_dh <- -u_t_u/((u_t_u+h)^2)
    gradient = c(dg_duks,dg_dh)

    fixed_covar = subthreshold_genes_h2$intercept*diag(nrow(sumstats_sig))
    S = matrix(data = 0, nrow = nrow(sumstats_sig)+1, ncol = nrow(sumstats_sig)+1)
    S[1:nrow(sumstats_sig), 1:nrow(sumstats_sig)] <- fixed_covar
    S[nrow(S), nrow(S)] <- subthreshold_genes_h2$heritabilities[2,ncol(merged_annotations)]^2
    S_out <- S
    gradient_out <- gradient
    frac_sig_se = sqrt(t(gradient) %*% S %*% gradient)
  }

  #If no genome-wide significant genes, then set mixed results to same as subthreshold genes model
  if (nrow(sumstats_sig) == 0){
    mixed_results <- subthreshold_genes_h2
    number_genomesig_genes <- 0
    genomesig_genes <- NA
    frac_sig = NA
    frac_sig_se = NA
  }

  #calculate the attenuation ratio

  a = subthreshold_genes_h2$intercept - (1/sumstats_true$N[1])
  mean_gammasq = mean(sumstats_true$gamma_sq)
  attenuation_ratio = a/mean_gammasq

  #standard error with delta method
  dg_da = 1/mean_gammasq
  dg_dgammasq = -(a - (1/sumstats_true$n_eff[1]))/(mean_gammasq^2)

  attenuation_gradient = matrix(data = c(dg_da,dg_dgammasq), nrow = 2, ncol = 1)


  jackknife_intercept = subthreshold_genes_h2$jackknife_betas[nrow(subthreshold_genes_h2$jackknife_betas),]
  jackknife_meangammasq= sapply(1:num_blocks, function(x) mean((sumstats_true$gamma_sq)[block_true != x]))
  jackknife_g_mat = cbind(jackknife_intercept - mean(jackknife_intercept), jackknife_meangammasq - mean(jackknife_meangammasq))
  S = ((num_blocks-1)/num_blocks)*(t(jackknife_g_mat) %*% jackknife_g_mat)

  attenuation_ratio_se = sqrt(t(attenuation_gradient) %*% S %*% attenuation_gradient)

  #calculate lambda gc
  chisq = (sumstats_true$N - 1) * (sumstats_true$gamma_sq)
  median_chisq = quantile(chisq, 0.5)
  lambda_gc = median_chisq/qchisq(0.5,1)

  jackknife_lambdagc = sapply(1:num_blocks, function(x) {
    quantile(chisq[block_true != x],0.5)/qchisq(0.5,1)
  })
  lambda_gc_se = sqrt(((num_blocks - 1)/num_blocks)*sum((jackknife_lambdagc - mean(jackknife_lambdagc))^2))

  if (all_models == FALSE){
    output <- list(mixed_model = mixed_results,
                   significant_genes = list(number_significant_genes = number_genomesig_genes,
                                            significant_genes = genomesig_genes,
                                            fraction_burdenh2_significant = frac_sig,
                                            fraction_burdenh2_significant_se = frac_sig_se),
                   qc = list(intercept = subthreshold_genes_h2$intercept,
                             intercept_se = subthreshold_genes_h2$intercept_se,
                             attenuation_ratio = attenuation_ratio,
                             attenuation_ratio_se = attenuation_ratio_se,
                             lambda_gc = lambda_gc,
                             lambda_gc_se = lambda_gc_se,
                             mu_genome = mu_genome))

  } else {
    output <- list(subthreshold_genes = subthreshold_genes_h2,
                   mixed_model = mixed_results,
                   significant_genes = list(number_significant_genes = number_genomesig_genes,
                                            significant_genes = genomesig_genes,
                                            fraction_burdenh2_significant = frac_sig,
                                            fraction_burdenh2_significant_se = frac_sig_se),
                   qc = list(intercept = subthreshold_genes_h2$intercept,
                             intercept_se = subthreshold_genes_h2$intercept_se,
                             attenuation_ratio = attenuation_ratio,
                             attenuation_ratio_se = attenuation_ratio_se,
                             lambda_gc = lambda_gc,
                             lambda_gc_se = lambda_gc_se,
                             mu_genome = mu_genome))
  }

  return(output)

}
