randomeffects_jackknife <- function(sumstats, merged_annotations,num_blocks,block,genomewide_correction,output_jackknife_h2, overdispersion, slope_correction, bivariate = FALSE, intercept = TRUE){
  #Estimate h2 using all genes
  if (overdispersion == TRUE){
    if (intercept == TRUE) {
      X = as.matrix(cbind(sumstats$burden_score*sumstats[colnames(merged_annotations[-1])],sumstats$burden_score,sumstats$overdispersion, rep(1,nrow(sumstats))))
      y = as.vector(sumstats$gamma_sq)
    } else {
      X = as.matrix(cbind(sumstats$burden_score*sumstats[colnames(merged_annotations[-1])],sumstats$burden_score,sumstats$overdispersion))
      y = as.vector(sumstats$gamma_sq)  - (1/sumstats$N)
    }
  } else {
    
    if (intercept == TRUE) {
      X = as.matrix(cbind(sumstats$burden_score*sumstats[colnames(merged_annotations[-1])],sumstats$burden_score, rep(1,nrow(sumstats))))
      y = as.vector(sumstats$gamma_sq)
    } else {
      X = as.matrix(cbind(sumstats$burden_score*sumstats[colnames(merged_annotations[-1])],sumstats$burden_score))
      y = as.vector(sumstats$gamma_sq) - (1/sumstats$N)
    }

  }

  #Compute per JK block "numerators" and "denominators" for calculation of OLS betas for all_block and leave-one-out
  jackknife_numerators = sapply(1:num_blocks, function(b,y,X) {return(t(X[block == b,])%*%y[block == b])}, y, X)
  jackknife_denominators = sapply(1:num_blocks, function(b,X) {return(t(X[block == b,])%*%X[block == b,])}, X)

  #Step 2.1: Compute total and annotation h2 from all_blocks
  all_block_beta_results <- solve(matrix(rowSums(jackknife_denominators),
                                         ncol = sqrt(length(rowSums(jackknife_denominators))),
                                         nrow = sqrt(length(rowSums(jackknife_denominators)))),
                                  rowSums(jackknife_numerators))
  if ((!is.null(slope_correction)) & (overdispersion == TRUE)){
    all_block_beta_results <- all_block_beta_results - c(rep(0,ncol(merged_annotations)-1), slope_correction, 0, 0)
  } else if ((!is.null(slope_correction)) & (overdispersion == FALSE)){
    all_block_beta_results <- all_block_beta_results - c(rep(0,ncol(merged_annotations)-1), slope_correction, 0)
  } else {
    all_block_beta_results <- all_block_beta_results
  }

  #now that the jackknife numerators and denominators are done, we can remove the null burden statistics

  sumstats_true = sumstats[sumstats$true,]
  block_true = block[sumstats$true]

  all_block_gene_h2_results <- as.matrix(cbind(sumstats_true[colnames(merged_annotations[-1])], 1)) %*% all_block_beta_results[1:ncol(merged_annotations)]*sumstats_true[,"burden_score"]
  h2 <- (t(as.matrix(cbind(sumstats_true[colnames(merged_annotations[-1])], 1))) %*% as.vector(all_block_gene_h2_results))
  #Step 2.2: Compute SEs using leave-1-out blocks
  jackknife_betas <- sapply(1:num_blocks,
                            function(b, numerators, denominators) {
                              lbo_numerators = numerators[,-b]
                              lbo_denominators = denominators[,-b]
                              solve(matrix(rowSums(lbo_denominators),
                                           ncol = sqrt(nrow(lbo_denominators)),
                                           nrow = sqrt(nrow(lbo_denominators))),
                                    rowSums(lbo_numerators))},
                            jackknife_numerators,
                            jackknife_denominators)

  if ((!is.null(slope_correction)) & (overdispersion == TRUE)){
    jackknife_betas <- sweep(jackknife_betas, MARGIN = 1, FUN = "-", STATS = c(rep(0,ncol(merged_annotations)-1), slope_correction, 0, 0))
  } else if ((!is.null(slope_correction)) & (overdispersion == FALSE)){
    jackknife_betas <- sweep(jackknife_betas, MARGIN = 1, FUN = "-", STATS = c(rep(0,ncol(merged_annotations)-1), slope_correction, 0))
  } else {
    jackknife_betas <- jackknife_betas
  }

  jackknife_geneh2 <- sapply(1:num_blocks,
                             function(b, jk_betas, sumstats_true, merged_annotations) {
                               (as.matrix(cbind(sumstats_true[colnames(merged_annotations[-1])], 1)[block_true != b,]) %*% jk_betas[1:ncol(merged_annotations),b]*sumstats_true[block_true != b,"burden_score"])
                             },
                             jackknife_betas,
                             sumstats_true,
                             merged_annotations, simplify = FALSE)
  if (intercept) {
    intercept_se = sqrt(((num_blocks - 1)/num_blocks)*sum((jackknife_betas[nrow(jackknife_betas),] - mean(jackknife_betas[nrow(jackknife_betas),]))^2))
  } else {
    intercept_se = NA
  }
  jackknife_annot_h2 <- sapply(1:num_blocks,
                               function(b, jk_geneh2s, sumstats_true, merged_annotations) {
                                 t(as.matrix(cbind(sumstats_true[colnames(merged_annotations[-1])], 1)[block_true != b,])) %*% as.vector(jk_geneh2s[[b]])
                               },
                               jackknife_geneh2,
                               sumstats_true,
                               merged_annotations)

  h2_se = sqrt(((num_blocks - 1)/num_blocks)*rowSums((jackknife_annot_h2 - rowMeans(jackknife_annot_h2))^2))
  #Step 2.3: Compute fraction of heritability explained by each gene set
  all_block_fractionh1_annot = sapply(1:(nrow(h2)-1), function(x) h2[x,1]/h2[nrow(h2),1])
  jackknife_fractionh2_annot <- t(sapply(1:(ncol(merged_annotations)-1),
                                         function(x) jackknife_annot_h2[x,]/jackknife_annot_h2[nrow(jackknife_annot_h2),]))
  fractionh2_annot_se = sqrt(((num_blocks - 1)/num_blocks)*rowSums((jackknife_fractionh2_annot - rowMeans(jackknife_fractionh2_annot))^2))

  h2_output_table <- rbind(t(h2),h2_se)
  colnames(h2_output_table)[ncol(h2_output_table)] = "total"
  rownames(h2_output_table) = c("h2","h2_se")

  fraction_output_table <- rbind(all_block_fractionh1_annot,fractionh2_annot_se)
  rownames(fraction_output_table) <- c("fraction_h2","fraction_h2_se")
  #Compute fraction of burden statistics per annotation
  fraction_burden_statistic <- (sumstats_true$burden_score %*% data.matrix(sumstats_true[colnames(merged_annotations[-1])]))/sum(sumstats_true$burden_score)
  enrichments <- data.frame(sweep(fraction_output_table, 2, fraction_burden_statistic, `/`))
  rownames(enrichments) <- c('enrichments', 'enrichment_se')

  if (bivariate){
    genomewide_h2 = ((sum(sumstats_true$w_t_beta)/sum(sumstats_true$burden_score)) - all_block_beta_results[length(all_block_beta_results)])
  } else{
    genomewide_h2 = ((sum(sumstats_true$w_t_beta)^2/sum(sumstats_true$burden_score)) - all_block_beta_results[length(all_block_beta_results)])

  }
  
  if (intercept) {
    intercept = all_block_beta_results[length(all_block_beta_results)]
  } else {
    intercept = NA
  }

  if (output_jackknife_h2) {
    return(list(heritabilities = h2_output_table,
                fractions = fraction_output_table,
                intercept = intercept,
                intercept_se = intercept_se,
                jackknife_betas = jackknife_betas,
                jackknife_h2 = jackknife_annot_h2,
                genomewide_h2 = genomewide_h2,
                enrichments_final = enrichments))
  } else {
    return(list(heritabilities = h2_output_table,
                fractions = fraction_output_table,
                intercept = intercept,
                intercept_se = intercept_se,
                genomewide_h2 = genomewide_h2,
                enrichments_final = enrichments))
  }
}
