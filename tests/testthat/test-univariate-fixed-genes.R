make_gene_sumstats <- function(n_genes = 120, n = 1000) {
  genes <- paste0("gene", seq_len(n_genes))
  fixed_gene_index <- 7
  gamma_sq <- rep(1e-5, n_genes)
  gamma_sq[fixed_gene_index] <- 0.05
  burden_score <- rep(0.01, n_genes)
  sample_size <- n + seq_len(n_genes)

  list(
    sumstats = data.frame(
      gene = genes,
      N = sample_size,
      gamma_sq = gamma_sq,
      w_t_beta = sqrt(gamma_sq * burden_score),
      burden_score = burden_score,
      overdispersion = rep(0.001, n_genes),
      chromosome = rep(1, n_genes),
      gene_position = seq_len(n_genes),
      num_variants = rep(1, n_genes)
    ),
    fixed_gene_index = fixed_gene_index,
    annotations = list(data.frame(
      gene = genes,
      even_gene = as.integer(seq_len(n_genes) %% 2 == 0)
    ))
  )
}

run_gene_level_bhr_h2 <- function(fixed_genes) {
  input <- make_gene_sumstats()

  suppressMessages(BHR_h2(
    sumstats = input$sumstats,
    annotations = input$annotations,
    num_blocks = 10,
    genomewide_correction = FALSE,
    fixed_genes = fixed_genes,
    output_jackknife_h2 = TRUE,
    overdispersion = FALSE,
    all_models = TRUE,
    num_null_conditions = 0,
    slope_correction = NULL,
    gwc_exclusion = NULL,
    intercept = FALSE,
    log = list(),
    start_time = Sys.time()
  ))
}

test_that("intercept = FALSE mixed model includes automatic fixed-gene h2", {
  input <- make_gene_sumstats()
  output <- run_gene_level_bhr_h2(fixed_genes = NULL)

  subthreshold_total <- output$subthreshold_genes$heritabilities["h2", "total"]
  subthreshold_total_se <- output$subthreshold_genes$heritabilities["h2_se", "total"]
  mixed_total <- output$mixed_model$heritabilities["h2", "total"]
  mixed_total_se <- output$mixed_model$heritabilities["h2_se", "total"]
  fixed_gene_index <- input$fixed_gene_index
  expected_fixed_h2 <- input$sumstats$gamma_sq[fixed_gene_index] - (1 / input$sumstats$N[fixed_gene_index])
  expected_fixed_var <- (4 * input$sumstats$gamma_sq[fixed_gene_index] / input$sumstats$N[fixed_gene_index]) +
    (4 / input$sumstats$N[fixed_gene_index]^2)

  expect_equal(output$significant_genes$number_significant_genes, 1)
  expect_equal(output$significant_genes$significant_genes, input$sumstats$gene[fixed_gene_index])
  expect_equal(
    as.numeric(mixed_total),
    as.numeric(subthreshold_total + expected_fixed_h2),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(mixed_total_se),
    as.numeric(sqrt(subthreshold_total_se^2 + expected_fixed_var)),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(output$significant_genes$fraction_burdenh2_significant),
    as.numeric(expected_fixed_h2 / mixed_total),
    tolerance = 1e-10
  )
  expect_equal(output$significant_genes$sig_table$varexplained, expected_fixed_h2)
  expect_equal(
    as.numeric(output$significant_genes$sig_table$fraction_burdenh2),
    as.numeric(expected_fixed_h2 / mixed_total),
    tolerance = 1e-10
  )
})

test_that("fixed_genes = character(0) adds no fixed-gene h2", {
  output <- run_gene_level_bhr_h2(fixed_genes = character(0))

  expect_equal(output$significant_genes$number_significant_genes, 0)
  expect_equal(
    as.numeric(output$mixed_model$heritabilities["h2", "total"]),
    as.numeric(output$subthreshold_genes$heritabilities["h2", "total"]),
    tolerance = 1e-10
  )
})
