library(bhr)
source("correct_winners_curse.R")
library(pracma)

########################################BHR basic estimates (h2, intercept) #############################

path = "/Users/daniel/Desktop/rare_h2/ms/genebass_summary_statistics_null/"
phenotype_file <- read.csv("~/rv_h2/reference_files/ms_phenotype_sheet.csv")
traits = phenotype_file[phenotype_file$phenotype_rg == 1,"phenotype_key"]
baseline_model <- read.table("~/rv_h2/reference_files/ms_baseline_oe5.txt")

bhr_summary_statistic_names <- c("bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0.0001_high0.001_group3",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_synonymous_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_synonymous_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_synonymous_low0.0001_high0.001_group3",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-benign_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-benign_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-benign_low0.0001_high0.001_group3",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0.0001_high0.001_group3")

results_holder <- matrix(data = NA, ncol = 23, nrow = length(bhr_summary_statistic_names) * length(traits))
counter = 1
for (ss in bhr_summary_statistic_names){
  summary_statistics <- readRDS(paste0(path,ss,".ms.munged.Rds"))
  for (trait in traits){
    n_bhr_trait <- head(summary_statistics[summary_statistics$phenotype_key == trait,"N"],1)
    print(c(trait, n_bhr_trait))
    output = BHR(mode = "univariate",
                 trait1_sumstats = summary_statistics[summary_statistics$phenotype_key == trait,],
                 annotations = list(baseline_model))
    
    results_holder[counter,] <- c(trait, 
                                  ss, 
                                  output$mixed_model$heritabilities[1,ncol(baseline_model)],
                                  output$mixed_model$heritabilities[2,ncol(baseline_model)],
                                  output$mixed_model$enrichments[1,1],
                                  output$mixed_model$enrichments[2,1],
                                  output$mixed_model$enrichments[1,2],
                                  output$mixed_model$enrichments[2,2],
                                  output$mixed_model$enrichments[1,3],
                                  output$mixed_model$enrichments[2,3],
                                  output$mixed_model$enrichments[1,4],
                                  output$mixed_model$enrichments[2,4],
                                  ((1-sum(output$mixed_model$fractions[1,]))/(1-sum(output$mixed_model$fraction_burden_score))),
                                  output$significant_genes$number_significant_genes,
                                  output$significant_genes$fraction_burdenh2_significant,
                                  output$significant_genes$fraction_burdenh2_significant_se,
                                  output$qc$intercept,
                                  output$qc$intercept_se,
                                  output$qc$attenuation_ratio,
                                  output$qc$attenuation_ratio_se,
                                  output$qc$lambda_gc,
                                  output$qc$lambda_gc_se,
                                  output$qc$mu_genome)
    
    print(paste0("Finished BHR estimate for ",trait," in summary statistic group ",ss)) 
    counter = counter + 1
  }
}
results_holder_df = as.data.frame(results_holder)
results_holder_df[,3:23] <- sapply(results_holder_df[,3:23],as.numeric)
colnames(results_holder_df) <- c("phenotype_key", "summary_statistic", "bhr_h2", "bhr_h2_se", 
                                 "bhr_enrichment_oe1", "bhr_enrichment_oe1_se", 
                                 "bhr_enrichment_oe2", "bhr_enrichment_oe2_se", 
                                 "bhr_enrichment_oe3", "bhr_enrichment_oe3_se", 
                                 "bhr_enrichment_oe4", "bhr_enrichment_oe4_se", 
                                 "bhr_enrichment_oe5",
                                 "n_significant_genes","fraction_h2_significant_genes", "fraction_h2_significant_genes_se",
                                 "intercept", "intercept_se", "attenuation_ratio", "attenuation_ratio_se", "lambda_gc", "lambda_gc_se", "mu_genome")
results_holder_df = merge(results_holder_df, phenotype_file[c("phenotype_key", "display_name", "phenocode", "n_bhr", "phenotype_core")], by = "phenotype_key")

write.csv(results_holder_df, "~/rv_h2/outputs/bhr_output_h2.csv")


######################Basic BHR estimates, plus slope correction ##########################
source("~/rv_h2/BHR.R")
source("~/rv_h2/BHR_h2.R")
source("~/rv_h2/randomeffects_jackknife.R")
source("~/rv_h2/BHR_meta.R")

path = "/Users/daniel/Desktop/rare_h2/ms/genebass_summary_statistics_null/"
phenotype_file <- read.csv("~/rv_h2/reference_files/ms_phenotype_sheet.csv")
traits = phenotype_file[phenotype_file$phenotype_rg == 1,"phenotype_key"]
baseline_model <- read.table("~/rv_h2/reference_files/ms_baseline_oe5.txt")

bhr_summary_statistic_names <- c("bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0.0001_high0.001_group3",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_synonymous_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_synonymous_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_synonymous_low0.0001_high0.001_group3",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-benign_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-benign_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-benign_low0.0001_high0.001_group3",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0.0001_high0.001_group3")

results_holder <- matrix(data = NA, ncol = 23, nrow = length(bhr_summary_statistic_names) * length(traits))
counter = 1
for (ss in bhr_summary_statistic_names){
  summary_statistics <- readRDS(paste0(path,ss,".ms.munged.Rds"))
  for (trait in traits){
    n_bhr_trait <- head(summary_statistics[summary_statistics$phenotype_key == trait,"N"],1)
    print(c(trait, n_bhr_trait))
    output = BHR(mode = "univariate",
                 trait1_sumstats = summary_statistics[summary_statistics$phenotype_key == trait,],
                 annotations = list(baseline_model),
                 slope_correction = 4.55087151/n_bhr_trait)
    
    results_holder[counter,] <- c(trait, 
                                  ss, 
                                  output$mixed_model$heritabilities[1,ncol(baseline_model)],
                                  output$mixed_model$heritabilities[2,ncol(baseline_model)],
                                  output$mixed_model$enrichments[1,1],
                                  output$mixed_model$enrichments[2,1],
                                  output$mixed_model$enrichments[1,2],
                                  output$mixed_model$enrichments[2,2],
                                  output$mixed_model$enrichments[1,3],
                                  output$mixed_model$enrichments[2,3],
                                  output$mixed_model$enrichments[1,4],
                                  output$mixed_model$enrichments[2,4],
                                  ((1-sum(output$mixed_model$fractions[1,]))/(1-sum(output$mixed_model$fraction_burden_score))),
                                  output$significant_genes$number_significant_genes,
                                  output$significant_genes$fraction_burdenh2_significant,
                                  output$significant_genes$fraction_burdenh2_significant_se,
                                  output$qc$intercept,
                                  output$qc$intercept_se,
                                  output$qc$attenuation_ratio,
                                  output$qc$attenuation_ratio_se,
                                  output$qc$lambda_gc,
                                  output$qc$lambda_gc_se,
                                  output$qc$mu_genome)
    
    print(paste0("Finished BHR estimate for ",trait," in summary statistic group ",ss)) 
    counter = counter + 1
  }
}
results_holder_df = as.data.frame(results_holder)
results_holder_df[,3:23] <- sapply(results_holder_df[,3:23],as.numeric)
colnames(results_holder_df) <- c("phenotype_key", "summary_statistic", "bhr_h2", "bhr_h2_se", 
                                 "bhr_enrichment_oe1", "bhr_enrichment_oe1_se", 
                                 "bhr_enrichment_oe2", "bhr_enrichment_oe2_se", 
                                 "bhr_enrichment_oe3", "bhr_enrichment_oe3_se", 
                                 "bhr_enrichment_oe4", "bhr_enrichment_oe4_se", 
                                 "bhr_enrichment_oe5",
                                 "n_significant_genes","fraction_h2_significant_genes", "fraction_h2_significant_genes_se",
                                 "intercept", "intercept_se", "attenuation_ratio", "attenuation_ratio_se", "lambda_gc", "lambda_gc_se", "mu_genome")
results_holder_df = merge(results_holder_df, phenotype_file[c("phenotype_key", "display_name", "phenocode", "n_bhr", "phenotype_core")], by = "phenotype_key")

write.csv(results_holder_df, "~/rv_h2/outputs/bhr_output_h2_ldcorrection.csv")

#############################################BHR gene set enrichments#####################
source("~/rv_h2/BHR.R")
source("~/rv_h2/BHR_h2.R")
source("~/rv_h2/randomeffects_jackknife.R")
source("~/rv_h2/BHR_meta.R")

phenotype_file <- read.csv("~/rv_h2/reference_files/ms_phenotype_sheet.csv")
gene_set_path = "/Users/daniel/Desktop/rare_h2/ms/gene_sets/gene_set_output/"
bhr_gene_sets <- c("brain_GABAergic",
                   "brain_Glutamatergic",
                   "cosmic_all",
                   "cosmic_oncogene",
                   "cosmic_tsg",
                   "ICA_cordblood_Erythroid",
                   "ICA_cordblood_Megakaryocytes",
                   "liver_Epithelial",
                   "segblood",
                   "segcortex",
                   "segliver",
                   "siggene_50NA",
                   "siggene_2453NA",
                   "siggene_3063NA",
                   "siggene_21001NA",
                   "siggene_30010NA",
                   "siggene_30080NA",
                   "siggene_30620NA",
                   "siggene_30680NA",
                   "siggene_30750NA",
                   "siggene_30770NA",
                   "siggene_30780NA",
                   "siggene_50NA")

summary_statistics <- readRDS("/Users/daniel/Desktop/rare_h2/ms/genebass_summary_statistics_null/bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1.ms.munged.Rds")
traits = phenotype_file[phenotype_file$phenotype_core == 1,"phenotype_key"]
baseline_model <- read.table("~/rv_h2/reference_files/ms_baseline_oe5.txt")

results_holder <- matrix(data = NA, ncol = 6, nrow = length(bhr_gene_sets) * length(traits))
counter = 1
for (bhr_gene_set in bhr_gene_sets){
  bhr_gene_set_annotation <- read.table(paste0(gene_set_path,"bhr_ms_" ,bhr_gene_set,".txt"), header = TRUE)
  for (trait in traits){
    output = BHR(mode = "univariate",
                 trait1_sumstats = summary_statistics[summary_statistics$phenotype_key == trait,],
                 annotations = list(baseline_model, bhr_gene_set_annotation))
    
    results_holder[counter,] <- c(trait, 
                                  bhr_gene_set, 
                                  output$mixed_model$fractions[1,ncol(baseline_model)],
                                  output$mixed_model$fractions[2,ncol(baseline_model)],
                                  output$mixed_model$enrichments[1,ncol(baseline_model)],
                                  output$mixed_model$enrichments[2,ncol(baseline_model)])
    counter = counter + 1
    print(paste0("Completed trait ",trait))
  }
  print(paste0("Completed gene set: ", bhr_gene_set))
}
results_holder_df = as.data.frame(results_holder)
results_holder_df[,3:6] <- sapply(results_holder_df[,3:6],as.numeric)
colnames(results_holder_df) <- c("phenotype_key", "gene_set", "fraction_h2", "fraction_h2_se", "enrichment", "enrichment_se")
results_holder_df$enrichment_z <- (results_holder_df$enrichment - 1) / results_holder_df$enrichment_se
results_holder_df = merge(results_holder_df, phenotype_file[c("phenotype_key", "display_name", "phenocode", "n_bhr", "phenotype_core")], by = "phenotype_key")
results_holder_df_rename <- results_holder_df
results_holder_df_rename$gene_set_display <- ifelse(results_holder_df_rename$gene_set == "brain_Glutamatergic", "Glutamatergic_neurons",
                                                    ifelse(results_holder_df_rename$gene_set == "segliver", "Liver",
                                                           ifelse(results_holder_df_rename$gene_set == "ICA_cordblood_Megakaryocytes", "Megakaryocytes",
                                                                  ifelse(results_holder_df_rename$gene_set == "cosmic_all", "Cancer_genes",
                                                                         ifelse(results_holder_df_rename$gene_set == "cosmic_oncogene", "Cancer_oncogenes",
                                                                                ifelse(results_holder_df_rename$gene_set == "segblood", "Whole_blood",
                                                                                       ifelse(results_holder_df_rename$gene_set == "cosmic_tsg", "Cancer_TSG",
                                                                                              ifelse(results_holder_df_rename$gene_set == "ICA_cordblood_Erythroid", "Erythrocyte",
                                                                                                     ifelse(results_holder_df_rename$gene_set == "brain_GABAergic", "GABAergic_neuron",
                                                                                                            ifelse(results_holder_df_rename$gene_set == "liver_Epithelial", "Hepatocyte",
                                                                                                                   ifelse(results_holder_df_rename$gene_set == "segcortex", "Cortex", NA)))))))))))

write.csv(results_holder_df_rename[c("gene_set", "gene_set_display", "phenotype_key", "display_name",
                                     "fraction_h2", "fraction_h2_se", "enrichment", "enrichment_se", "enrichment_z")], "~/rv_h2/outputs/bhr_output_geneset.csv")

#################################################BHR aggregate for each trait, and for all traits############################################
source("~/rv_h2/BHR.R")
source("~/rv_h2/BHR_h2.R")
source("~/rv_h2/randomeffects_jackknife.R")
source("~/rv_h2/BHR_meta.R")

path = "/Users/daniel/Desktop/rare_h2/ms/genebass_summary_statistics_null/"
phenotype_file <- read.csv("~/rv_h2/reference_files/ms_phenotype_sheet.csv")
traits = phenotype_file[phenotype_file$phenotype_rg == 1,"phenotype_key"]
baseline_model <- read.table("~/rv_h2/reference_files/ms_baseline_oe5.txt")

bhr_summary_statistic_names <- c("bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0.0001_high0.001_group3",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0.0001_high0.001_group3")

#Aggregate h2 for each trait for g1-g3 (ultra-rare + rare)
ss_list <- list(readRDS(paste0(path, bhr_summary_statistic_names[1],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[2],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[3],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[4],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[5],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[6],".ms.munged.Rds")))

results_holder <- matrix(data = NA, ncol = 3, nrow = length(traits))
counter = 1
for (trait in traits){
  output <- BHR(mode = 'aggregate', 
                ss_list_trait1 = ss_list, trait_list = list(trait), 
                annotations = baseline_model)
  results_holder[counter,] <- c(trait,
                                output$aggregated_mixed_model_h2,
                                output$aggregated_mixed_model_h2se)
  print(paste0("Completed ",trait))
  counter = counter + 1
}
results_holder_df = as.data.frame(results_holder)
results_holder_df[,2:3] <- sapply(results_holder_df[,2:3],as.numeric)
colnames(results_holder_df) <- c("phenotype_key", "aggregated_h2", "aggregate_h2_se")
results_holder_df = merge(results_holder_df, phenotype_file[c("phenotype_key", "display_name")], by = "phenotype_key")
results_holder_df$ss_group = "g13_lof_mis"

write.csv(results_holder_df, "~/rv_h2/outputs/bhr_output_aggregated_per_trait_h2.csv")


#Aggregate across all core traits
results_holder <- matrix(data = NA, ncol = 3, nrow = 1)
output <- BHR(mode = 'aggregate', ss_list_trait1 = ss_list, 
              trait_list = phenotype_file[phenotype_file$phenotype_core == 1,"phenotype_key"], 
              annotations = baseline_model)
results_holder[1,] <- c("All_traits",
                        output$aggregated_mixed_model_h2,
                        output$aggregated_mixed_model_h2se)
results_holder_df = as.data.frame(results_holder)
results_holder_df[,2:3] <- sapply(results_holder_df[,2:3],as.numeric)
colnames(results_holder_df) <- c("phenotype_key", "aggregated_h2", "aggregated_h2_se")
write.csv(results_holder_df, "~/rv_h2/outputs/bhr_output_aggregated_all_traits_h2.csv")



#################################################BHR aggregate for each trait, and for all traits, with slope correction ############################################
source("~/rv_h2/BHR.R")
source("~/rv_h2/BHR_h2.R")
source("~/rv_h2/randomeffects_jackknife.R")
source("~/rv_h2/BHR_meta.R")

path = "/Users/daniel/Desktop/rare_h2/ms/genebass_summary_statistics_null/"
phenotype_file <- read.csv("~/rv_h2/reference_files/ms_phenotype_sheet.csv")
traits = phenotype_file[phenotype_file$phenotype_rg == 1,"phenotype_key"]
baseline_model <- read.table("~/rv_h2/reference_files/ms_baseline_oe5.txt")

bhr_summary_statistic_names <- c("bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0.0001_high0.001_group3",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0.0001_high0.001_group3")

#Aggregate h2 for each trait for g1-g3 (ultra-rare + rare)
ss_list <- list(readRDS(paste0(path, bhr_summary_statistic_names[1],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[2],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[3],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[4],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[5],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[6],".ms.munged.Rds")))

results_holder <- matrix(data = NA, ncol = 3, nrow = length(traits))
counter = 1
summary_statistics <- readRDS(paste0(path,"bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1.ms.munged.Rds"))
for (trait in traits){
  n_bhr_trait <- head(summary_statistics[summary_statistics$phenotype_key == trait,"N"],1)
  output <- BHR(mode = 'aggregate', 
                ss_list_trait1 = ss_list, trait_list = list(trait), 
                annotations = baseline_model,
                slope_correction = 4.55087151/n_bhr_trait)
  results_holder[counter,] <- c(trait,
                                output$aggregated_mixed_model_h2,
                                output$aggregated_mixed_model_h2se)
  print(paste0("Completed ",trait))
  counter = counter + 1
}
results_holder_df = as.data.frame(results_holder)
results_holder_df[,2:3] <- sapply(results_holder_df[,2:3],as.numeric)
colnames(results_holder_df) <- c("phenotype_key", "aggregated_h2", "aggregate_h2_se")
results_holder_df = merge(results_holder_df, phenotype_file[c("phenotype_key", "display_name")], by = "phenotype_key")
results_holder_df$ss_group = "g13_lof_mis"

write.csv(results_holder_df, "~/rv_h2/outputs/bhr_output_aggregated_per_trait_h2_ldcorrection.csv")


#Aggregate across all core traits
source("~/rv_h2/BHR.R")
source("~/rv_h2/BHR_h2.R")
source("~/rv_h2/randomeffects_jackknife.R")
source("~/rv_h2/BHR_meta.R")

path = "/Users/daniel/Desktop/rare_h2/ms/genebass_summary_statistics_null/"
phenotype_file <- read.csv("~/rv_h2/reference_files/ms_phenotype_sheet.csv")
baseline_model <- read.table("~/rv_h2/reference_files/ms_baseline_oe5.txt")

bhr_summary_statistic_names <- c("bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0.0001_high0.001_group3",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0_high1e-05_group1",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low1e-05_high0.0001_group2",
                                 "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0.0001_high0.001_group3")

#Aggregate h2 for each trait for g1-g3 (ultra-rare + rare)
ss_list <- list(readRDS(paste0(path, bhr_summary_statistic_names[1],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[2],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[3],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[4],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[5],".ms.munged.Rds")),
                readRDS(paste0(path, bhr_summary_statistic_names[6],".ms.munged.Rds")))
#remember to add slope line to meta for this
results_holder <- matrix(data = NA, ncol = 3, nrow = 1)
output <- BHR(mode = 'aggregate', ss_list_trait1 = ss_list, 
              trait_list = phenotype_file[phenotype_file$phenotype_core == 1,"phenotype_key"], 
              annotations = baseline_model)
results_holder[1,] <- c("All_traits",
                        output$aggregated_mixed_model_h2,
                        output$aggregated_mixed_model_h2se)
results_holder_df = as.data.frame(results_holder)
results_holder_df[,2:3] <- sapply(results_holder_df[,2:3],as.numeric)
colnames(results_holder_df) <- c("phenotype_key", "aggregated_h2", "aggregated_h2_se")
write.csv(results_holder_df, "~/rv_h2/outputs/bhr_output_aggregated_all_traits_h2_ldcorrection.csv")





####Significant gene fractions

phenotype_file <- read.csv("~/rv_h2/reference_files/ms_phenotype_sheet.csv")
sig_genes = read.csv("~/rv_h2/reference_files/siggenes_for_AMM.csv")
traits = phenotype_file[phenotype_file$phenotype_core == 1,"phenotype_key"]
summary_statistics_plof_grp1 <- readRDS("~/Documents/oconnor_rotation/rarevariantproject/bhr_ms_gene_ss_400k_withnullburden_pLoF_nvar449780_low0_high1e-05_group1.ms.munged.Rds")
baseline_model <- read.table("~/rv_h2/reference_files/ms_baseline_oe5.txt")
consensus_genes <- read.table("~/rv_h2/reference_files/bhr_ms_consensus_gene_list.txt")[c("gene_id", "gene")]

get_frac_assns <- function(trait){
  print(trait)
  trait_sumstats = summary_statistics_plof_grp1[summary_statistics_plof_grp1$phenotype_key == trait,]
  
  sig_genes_trait = sig_genes$gene[sig_genes$id == trait]
  print(length(sig_genes_trait))
  output = BHR(mode = "univariate",
               trait1_sumstats = trait_sumstats,
               annotations = list(baseline_model),
               fixed_genes = sig_genes_trait,
               gwc_exclusion = FALSE)
  
  trait_sumstats_sig = trait_sumstats[trait_sumstats$gene %in% sig_genes_trait,]
  
  trait_sig_df <- data.frame(phenotype_key = trait,
                             gene = trait_sumstats_sig$gene,
                             varexplained = trait_sumstats_sig$w_t_beta^2/trait_sumstats_sig$burden_score,
                             bhr_h2 = output$mixed_model$heritabilities[1,ncol(output$mixed_model$heritabilities)],
                             frac_sig = output$significant_genes$fraction_burdenh2_significant,
                             frac_sig_se = output$significant_genes$fraction_burdenh2_significant_se)
  trait_sig_df$chisq = trait_sumstats$N[1]*trait_sig_df$varexplained
  
  thresh = qchisq(p = 0.05/nrow(trait_sumstats),df = 1,lower.tail = FALSE)
  trait_sig_df$chisq_winnerscursecorr = correct_winners_curse(trait_sig_df$chisq,thresh)
  trait_sig_df$varexplained_winnerscursecorr = trait_sig_df$chisq_winnerscursecorr/trait_sumstats$N[1]
  trait_sig_df$frac_sig_winnerscursecorr = sum(trait_sig_df$varexplained_winnerscursecorr)/trait_sig_df$bhr_h2[1]
  return(trait_sig_df)
}

sig_results <- lapply(traits[traits %in% sig_genes$id],get_frac_assns)

sig_df = sig_results[[1]]
for (trait in 2:length(sig_results)){
  print(trait)
  sig_df = rbind(sig_df,sig_results[[trait]])
}

sig_df$display_name = phenotype_file$display_name[match(sig_df$phenotype_key,phenotype_file$phenotype_key)]
sig_df$proportion_bhr_explained = sig_df$varexplained_winnerscursecorr/sig_df$bhr_h2
sig_df$labelgene = consensus_genes$gene[match(sig_df$gene,consensus_genes$gene_id)]

ordered_sigdf = sig_df[order(sig_df$phenotype_key,-sig_df$varexplained_winnerscursecorr),]
write.csv(ordered_sigdf, "~/rv_h2/outputs/BHR_Significant_Genes_Info.csv", quote = FALSE, row.names = FALSE)


####Significant gene fractions, common variant space

variant_table <- fread2("~/Documents/oconnor_rotation/rarevariantproject/variants_MAFfilt_subsetcols.tsv", header = TRUE)
sumstat_lookup <- read.csv("reference_files/sumstat_lookup.csv", header = FALSE)

sigclump_df = data.frame()
for (trait in 1:nrow(sumstat_lookup)){
  print(trait)
  print(sumstat_lookup$V21[trait])
  
  print("loading files")
  clumpfile = paste0("~/Documents/oconnor_rotation/rarevariantproject/final_commonvar_manuscript/significant_clump/clump_output/", strsplit(sumstat_lookup$V21[trait],split = ".bgz|.gz")[[1]],".plink.tsv_clump.clumped")
  gwas_file = paste0("~/Documents/oconnor_rotation/rarevariantproject/final_commonvar_manuscript/significant_clump/gwas_sumstats/",strsplit(sumstat_lookup$V10[trait],split = ".bgz|.gz")[[1]],"_lean")
  ldsc_file = paste0("~/Documents/oconnor_rotation/rarevariantproject/final_commonvar_manuscript/significant_clump/input_sumstats/",strsplit(sumstat_lookup$V7[trait],split = ".bgz|.gz")[[1]])
  
  
  clumps = read.table(clumpfile, header = TRUE)
  clumps = clumps[clumps$P < 5e-8,]  
  
  gwas = fread2(gwas_file, header = FALSE)
  ldsc = read.table(ldsc_file, header = TRUE)
  names(gwas) <- c("variant","minor_AF","beta")
  
  print("extracting significant clumps variances")
  snps = clumps$SNP
  chrs = variant_table$chr[match(snps,variant_table$rsid)]
  bps = variant_table$pos[match(snps,variant_table$rsid)]
  refs = variant_table$ref[match(snps,variant_table$rsid)]
  alts = variant_table$alt[match(snps,variant_table$rsid)]
  
  identifiers = paste(chrs,bps,refs,alts,sep = ":")
  
  betas = gwas$beta[match(identifiers,gwas$variant)]
  mafs = gwas$minor_AF[match(identifiers,gwas$variant)]
  varexplained = (betas*(sqrt(2*mafs*(1 - mafs))))^2
  
  traitclumpdf <- data.frame(trait = sumstat_lookup$V21[trait],
                             SNP = snps,
                             varexplained = varexplained,
                             N = ldsc$N[1])
  traitclumpdf = traitclumpdf[!is.na(traitclumpdf$varexplained),]
  traitclumpdf$chisq = traitclumpdf$varexplained * ldsc$N[1]
  
  thresh = qchisq(5e-8,1,lower.tail = FALSE)
  
  traitclumpdf$chisq_winnerscursecorr = correct_winners_curse(traitclumpdf$chisq,thresh)
  traitclumpdf$varexplained_winnerscursecorr = traitclumpdf$chisq_winnerscursecorr/ldsc$N[1]
  
  sigclump_df <- rbind(sigclump_df,traitclumpdf)
  
}

#aside: get Ns

get_N <- function(trait){
  print(trait)
  print(sumstat_lookup$V21[trait])
  
  ldsc_file = paste0("~/Documents/oconnor_rotation/rarevariantproject/final_commonvar_manuscript/significant_clump/input_sumstats/",strsplit(sumstat_lookup$V7[trait],split = ".bgz|.gz")[[1]])
  
  ldsc = read.table(ldsc_file, header = TRUE)
  
  return(ldsc$N[1])
  
}
ldsc_h2 = read.csv("~/Documents/oconnor_rotation/rarevariantproject/rv_h2/ldsc_h2.csv")
ldsc_h2$N = sapply(1:nrow(ldsc_h2), get_N)
write.csv(ldsc_h2, "outputs//ldsc_h2.csv", row.names = FALSE, quote = FALSE)
sigclump_df <- sigclump_df[sigclump_df$trait %in% phenotype_file$phenotype_key[phenotype_file$phenotype_core ==1],]
sigclump_df <- sigclump_df[!is.na(sigclump_df$varexplained),]
get_HESSh2 <- function(trait){
  print(trait)
  HESS_h2_result = read.table(paste0("~/Documents/oconnor_rotation/rarevariantproject/final_commonvar_manuscript/HESS/step2_output/",  trait,"_HESSformat.tsv_step2.log_HESSh2"), sep = " ")
  
  return(as.numeric(HESS_h2_result$V5))
}

get_HESSh2_se <- function(trait){
  HESS_h2_result = read.table(paste0("~/Documents/oconnor_rotation/rarevariantproject/final_commonvar_manuscript/HESS/step2_output/",  trait,"_HESSformat.tsv_step2.log_HESSh2"), sep = " ")
  
  return(parse_number(HESS_h2_result$V6))
}

sigclump_df$HESS_h2 = sapply(sigclump_df$trait,get_HESSh2)
sigclump_df$HESS_h2_se = sapply(sigclump_df$trait,get_HESSh2_se)
sigclump_df$fraction_HESS = sigclump_df$varexplained_winnerscursecorr/sigclump_df$HESS_h2

sigclump_df_ordered = sigclump_df[order(sigclump_df$trait,-sigclump_df$varexplained_winnerscursecorr),]
names(sigclump_df_ordered)[1] <- "phenotype_key"

sigclump_df_ordered$fraction_HESS_sig = sapply(sigclump_df_ordered$phenotype_key,
                                               function(x) {
                                                 return(sum(sigclump_df_ordered$varexplained_winnerscursecorr[sigclump_df_ordered$phenotype_key == x])/(sigclump_df_ordered$HESS_h2[sigclump_df_ordered$phenotype_key == x][1]))
                                               })

write.csv(sigclump_df_ordered,"outputs/common_significant_associations.csv", quote = FALSE, row.names = FALSE)

####Significant gene fractions, common variant space (HESS partitions)
HESS_df <- data.frame()

for (trait in phenotype_file$phenotype_key[phenotype_file$phenotype_core == 1]){
  #print(trait)
  trait_hess_results <- read.table(paste0("~/Documents/oconnor_rotation/rarevariantproject/final_commonvar_manuscript/HESS/step2_output/",
                                          trait,
                                          "_HESSformat.tsv_step2.txt"),
                                   header = TRUE)
  
  trait_hess_results$trait = trait
  HESS_df = rbind(HESS_df,trait_hess_results)
  
}

write.csv(HESS_df,"outputs/HESS_output.csv", row.names = FALSE, quote = FALSE)
######################################## Genetic Correlation #############################
#Read in results from cluster

ldsc_rg_guide <- read.table("~/Documents/oconnor_rotation/rarevariantproject/final_commonvar_manuscript/LDSC_rg/ldsc_rg_guide.tsv", sep = ":")

get_ldsc_rg <- function(pair){
  input = read.table(paste0("~/Documents/oconnor_rotation/rarevariantproject/final_commonvar_manuscript/LDSC_rg/output/",pair,".txt"), sep = " ")
  return(as.numeric(input$V3[1]))
}

get_ldsc_rg_se <- function(pair){
  input = read.table(paste0("~/Documents/oconnor_rotation/rarevariantproject/final_commonvar_manuscript/LDSC_rg/output/",pair,".txt"), sep = " ")
  return(parse_number(input$V4[1]))
}

ldsc_rg_guide$ldsc_rg = sapply(1:nrow(ldsc_rg_guide), get_ldsc_rg)
ldsc_rg_guide$ldsc_rg_se = sapply(1:nrow(ldsc_rg_guide), get_ldsc_rg_se)


get_bhr_rg <- function(pair){
  print(pair)
  input = read.table(paste0("~/Documents/oconnor_rotation/rarevariantproject/BHR_rg_cluster/output/",pair,".txt"), header = TRUE)
  return(input$x[1])
}

get_bhr_rg_se <- function(pair){
  print(pair)
  input = read.table(paste0("~/Documents/oconnor_rotation/rarevariantproject/BHR_rg_cluster/output/",pair,".txt"), header = TRUE)
  return(input$x[2])
}

ldsc_rg_guide$bhr_rg = sapply(1:nrow(ldsc_rg_guide), get_bhr_rg)
ldsc_rg_guide$bhr_rg_se = sapply(1:nrow(ldsc_rg_guide), get_bhr_rg_se)

bhr_h2 = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "bhr_h2"))
bhr_h2_plofgrp1 =bhr_h2[bhr_h2$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1",]

ldsc_rg_guide$bhr_h2_trait1 = bhr_h2_plofgrp1$bhr_h2[match(ldsc_rg_guide$V3,bhr_h2_plofgrp1$phenotype_key)]
ldsc_rg_guide$bhr_h2_trait1_se = bhr_h2_plofgrp1$bhr_h2_se[match(ldsc_rg_guide$V3,bhr_h2_plofgrp1$phenotype_key)]
ldsc_rg_guide$bhr_h2_trait2 = bhr_h2_plofgrp1$bhr_h2[match(ldsc_rg_guide$V6,bhr_h2_plofgrp1$phenotype_key)]
ldsc_rg_guide$bhr_h2_trait2_se = bhr_h2_plofgrp1$bhr_h2_se[match(ldsc_rg_guide$V6,bhr_h2_plofgrp1$phenotype_key)]

ldsc_rg_guide$trait1_bhrh2_z = ldsc_rg_guide$bhr_h2_trait1/ldsc_rg_guide$bhr_h2_trait1_se 
ldsc_rg_guide$trait2_bhrh2_z = ldsc_rg_guide$bhr_h2_trait2/ldsc_rg_guide$bhr_h2_trait2_se 

rg_output = ldsc_rg_guide[,c("V3","V6","ldsc_rg","ldsc_rg_se","bhr_rg","bhr_rg_se","trait1_bhrh2_z","trait2_bhrh2_z")]
rg_output_sigbhrh2 = rg_output[abs(rg_output$trait1_bhrh2_z) > 1.96 & abs(rg_output$trait2_bhrh2_z) > 1.96,]
names(rg_output_sigbhrh2) <- c("trait1","trait2","ldsc_rg","ldsc_rg_se","bhr_rg","bhr_rg_se")
names(rg_output) <- c("trait1","trait2","ldsc_rg","ldsc_rg_se","bhr_rg","bhr_rg_se")


write.csv(rg_output,"outputs/rg_output.csv", quote = FALSE, row.names = FALSE)

#rg between missense and plof

plof_sumstats = readRDS("~/Documents/oconnor_rotation/rarevariantproject/bhr_ms_gene_ss_400k_withnullburden_pLoF_nvar449780_low0_high1e-05_group1.ms.munged.Rds")
missense_sumstats = readRDS("~/Documents/oconnor_rotation/rarevariantproject/bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_nvar1590551_low0_high1e-05_group1.ms.munged.Rds")

traits = phenotype_file$phenotype_key[phenotype_file$phenotype_core == 1]

get_plof_missense_stats <- function(trait){
  plof_sumstats_trait = plof_sumstats[plof_sumstats$phenotype_key == trait,]
  missense_sumstats_trait = missense_sumstats[missense_sumstats$phenotype_key == trait,]
  
  bhr_plof = BHR_h2(plof_sumstats_trait,
                    annotations = list(baseline_model), 
                    num_blocks = 100,
                    genomewide_correction = FALSE,
                    gwc_exclusion = NULL,
                    overdispersion = FALSE, 
                    num_null_conditions = 0,
                    output_jackknife_h2 = FALSE, 
                    fixed_genes = NULL,
                    all_models = FALSE,
                    slope_correction = FALSE )
  
  bhr_missense = BHR_h2(missense_sumstats_trait,
                        annotations = list(baseline_model), 
                        num_blocks = 100,
                        genomewide_correction = FALSE,
                        gwc_exclusion = NULL,
                        overdispersion = FALSE, 
                        num_null_conditions = 0,
                        output_jackknife_h2 = FALSE, 
                        fixed_genes = NULL,
                        all_models = FALSE,
                        slope_correction = FALSE )
  if (bhr_plof$mixed_model$heritabilities[1,5] < 0 | bhr_missense$mixed_model$heritabilities[1,5] < 0 ){
    output = list(plof_h2 = bhr_plof$mixed_model$heritabilities[1,5],
                  plof_h2_se = bhr_plof$mixed_model$heritabilities[2,5],
                  missense_h2 = bhr_missense$mixed_model$heritabilities[1,5],
                  missense_h2_se = bhr_missense$mixed_model$heritabilities[2,5],
                  rg = NA,
                  rg_se = NA)
    print(output)
    return(output)
  }
  
  bhr_rg_trait =  BHR_rg(plof_sumstats_trait, 
                         missense_sumstats_trait, 
                         annotations = list(baseline_model), 
                         num_blocks = 100,
                         genomewide_correction = FALSE,
                         overdispersion = FALSE, 
                         num_null_conditions = 0,
                         output_jackknife_rg = FALSE, 
                         fixed_genes = NULL)
  output = list(plof_h2 = bhr_plof$mixed_model$heritabilities[1,5],
                plof_h2_se = bhr_plof$mixed_model$heritabilities[2,5],
                missense_h2 = bhr_missense$mixed_model$heritabilities[1,5],
                missense_h2_se = bhr_missense$mixed_model$heritabilities[2,5],
                rg = bhr_rg_trait$rg$rg_mixed,
                rg_se = bhr_rg_trait$rg$rg_mixed_se)
  print(output)
  return(output)
}

plof_missense_compare_bhr <- sapply(traits, get_plof_missense_stats)
names = rownames(plof_missense_compare_bhr)
plof_missense_compare_bhr_df = as.data.frame(t(plof_missense_compare_bhr))
plof_missense_compare_bhr_df = as.data.frame(sapply(1:ncol(plof_missense_compare_bhr_df), function(x) as.numeric(unlist(plof_missense_compare_bhr_df[,x]))))
names(plof_missense_compare_bhr_df) <- names
plof_missense_compare_bhr_df$display_name = phenotype_file$display_name[phenotype_file$phenotype_core ==1]

write.csv(plof_missense_compare_bhr_df,"outputs/plof_missense_compare.csv", quote = FALSE, row.names = FALSE)


######################################## SCZ and BP #############################

#BipEx
#Read in the baseline model file
library(pracma)
library(tidyverse)
library(readr)
baseline_model <- read.table("~/Documents/oconnor_rotation/rarevariantproject/final_manuscript_repo/bhr/reference_files/ms_baseline_oe5.txt")

#Read in the publicly available bipex variant table.
bp_variantlevel_bipex <- bigreadr::fread2("~/Downloads/BipEx_variant_results.tsv")


#Subset variants 

variant_filter = bp_variantlevel_bipex$group == "Bipolar Disorder" & #Filter to Bipolar Disorder counts
  bp_variantlevel_bipex$in_analysis == TRUE & #Use variant filter from Palmer et al, 2022
  !is.na(bp_variantlevel_bipex$gene_id) & #Remove variants with NA gene ID 
  str_detect(bp_variantlevel_bipex$locus, "^chr\\d") #Subset to autosomal variants

bp_variantlevel_bipex <- bp_variantlevel_bipex[variant_filter,
                                               c("gene_id",
                                                 "consequence",
                                                 "ac_case",
                                                 "ac_ctrl",
                                                 "locus",
                                                 "mpc")]  
#Function to wrangle into BHR sumstats format
wrangle_sumstats <- function(table,n_cases,n_controls, var_filter) {
  
  #Filter to variants of interest
  table = table[var_filter,] 
  
  #Compute sample prevalence, will be used to compute per-sd beta
  prevalence = n_cases/(n_cases+n_controls) 
  
  #Compute variant MAF in cases, will be used to compute per-sd beta
  table$AF_case = table$ac_case/(2*n_cases) 
  
  #Compute variant MAFoverall, will be used to compute per-sd beta
  table$AF = (table$ac_case + table$ac_ctrl)/(2*(n_cases + n_controls)) 
  
  #calculate per-sd betas
  table$beta = (2 * (table$AF_case - table$AF) * prevalence)/sqrt(2 * table$AF * (1 - table$AF) * prevalence * (1 - prevalence))  
  
  #calculate variant variances
  table$twopq = 2*table$AF * (1 - table$AF) 
  
  #convert betas from per-sd (i.e. sqrt(variance explained)) to per-allele (i.e. in units of phenotype) betas.
  #per-allele are the usual betas reported by an exome wide association study.
  table$beta_perallele = table$beta/sqrt(table$twopq)
  
  #aggregate into gene-level table.
  sumstats = data.frame(gene = unique(table$gene_id))
  
  #each element of variant_variances is a list of variances for all variants in a gene.
  sumstats$variant_variances = lapply(sumstats$gene,
                                      function(x) table$twopq[table$gene_id == x]) 
  
  #each element of variant_variances is a list of per-allele betas for all variants in a gene.
  sumstats$betas = lapply(sumstats$gene,
                          function(x) table$beta_perallele[table$gene_id == x]) 
  
  
  names(sumstats) <- c("gene", "variant_variances","betas")
  
  #add chromosome and position information
  #position doesn't need to be super precise, as it is only used to order genes for jackknife
  sumstats$gene_position <- parse_number(sapply(strsplit(table$locus[match(sumstats$gene,table$gene_id)], split = ":"), function(x) x[[2]]))
  sumstats$chromosome = parse_number(sapply(strsplit(table$locus[match(sumstats$gene,table$gene_id)], split = ":"), function(x) x[[1]]))
  
  #N = sum of case and control counts
  sumstats$N = n_cases + n_controls
  
  #we have found that in these smaller sample analyses, there are some genes with
  #large burden scores that are clearly outliers
  #we remove genes with burden scores more than 8 sd from the mean as a conservative filter.
  burdenscores = sapply(sumstats$variant_variances, function(x) sum(x))
  sumstats <- sumstats[abs(scale(burdenscores)) < 8,]
  
  
  return(sumstats)
}

bp_sumstats_ptv = wrangle_sumstats(bp_variantlevel_bipex,
                                   14210, #N from Palmer et al, 2022
                                   14422, #N from Palmer et al, 2022
                                   bp_variantlevel_bipex$consequence == "ptv")

bp_sumstats_missenseMPC2 = wrangle_sumstats(bp_variantlevel_bipex[!is.na(bp_variantlevel_bipex$mpc),],
                                            14210,
                                            14422,
                                            bp_variantlevel_bipex$consequence[!is.na(bp_variantlevel_bipex$mpc)] %in% c("damaging_missense", "other_missense") &
                                              bp_variantlevel_bipex$mpc[!is.na(bp_variantlevel_bipex$mpc)] > 2)

bp_sumstats_synonymous = wrangle_sumstats(bp_variantlevel_bipex,
                                          14210,
                                          14422,
                                          bp_variantlevel_bipex$consequence == "synonymous")

#Run BHR
bp_ptv_bhr <- bhr::BHR(bp_sumstats_ptv, 
                       annotations = list(baseline_model), #baseline model including constraint annotations
                       num_blocks = 100, #number of blocks for jackknife
                       num_null_conditions = 5, #5*num_genes null moment conditions
                       mode = "univariate") #run in univariate mode to compute burden h2

bp_missense_bhr <- bhr::BHR(bp_sumstats_missenseMPC2, 
                                annotations = list(baseline_model),
                                num_blocks = 100,
                                num_null_conditions = 5,
                                mode = "univariate")

bp_synonymous_bhr <- bhr::BHR(bp_sumstats_synonymous, 
                              annotations = list(baseline_model),
                              num_blocks = 100,
                              genomewide_correction = FALSE, 
                              num_null_conditions = 5,
                              mode = "univariate")

#convert observed scale to liability scale h2
#calculate sample prevalence
bp_sample_prevalence = 14210/(14210 + 14422) 

#this population prevalence estimate is from Ferrari et al, 2016 Bipolar Disorders
population_prevalence_bp = 0.007 


obs2lia_factor <- function(K, P){
  X <- qnorm(K,lower.tail=FALSE)
  z <- (1/sqrt(2*pi))*(exp(-(X**2)/2))
  factor <- (K*(1-K)*K*(1-K))/(P*(1-P)*(z**2))
  return(factor)
}

bp_scalingfactor = obs2lia_factor(population_prevalence_bp,bp_sample_prevalence)

#Burden heritability of schizophrenia
#Read in table from SCHEMA website
SCHEMA_variants <- bigreadr::fread2("~/Downloads/SCHEMA_variant_results.tsv",
                                select = c("gene_id",
                                           "consequence",
                                           "ac_case",
                                           "ac_ctrl",
                                           "group",
                                           "locus",
                                           "in_analysis"))

#Filter variants
SCHEMA_variant_filter = (SCHEMA_variants$in_analysis |  SCHEMA_variants$consequence == "synonymous_variant") & #Either filtered SCHEMA damaging variant, or synonymous
  SCHEMA_variants$group %in% c("EUR (exomes)","EUR (gnomAD exomes)","EUR-N (exomes)") & #In one of the primary EUR cohorts
  !is.na(SCHEMA_variants$gene_id) &  #non-NA gene ID
  str_detect(SCHEMA_variants$locus, "X", negate = TRUE) & #autosomal
  str_detect(SCHEMA_variants$locus, "Y", negate = TRUE) #autosomal

SCHEMA_variants = SCHEMA_variants[SCHEMA_variant_filter,]
SCHEMA_variants = SCHEMA_variants[!is.na(SCHEMA_variants$ac_case),]

#Aggregate variants into nextera and non-nextera variant tables
EUR_Nextera <- SCHEMA_variants[SCHEMA_variants$group %in% c("EUR (exomes)","EUR (gnomAD exomes)"),]
EUR_agg_Nextera <- aggregate(cbind(EUR_Nextera$ac_case,EUR_Nextera$ac_ctrl), by = list(EUR_Nextera$locus),FUN = sum)
EUR_agg_Nextera <- EUR_agg_Nextera[EUR_agg_Nextera$V1 + EUR_agg_Nextera$V2 <= 5,]
EUR_agg_Nextera$gene_id = EUR_Nextera$gene_id[match(EUR_agg_Nextera$Group.1,EUR_Nextera$locus)]
EUR_agg_Nextera$consequence = EUR_Nextera$consequence[match(EUR_agg_Nextera$Group.1,EUR_Nextera$locus)]
names(EUR_agg_Nextera) <- c("locus","ac_case","ac_ctrl","gene_id","consequence")

EUR_NonNextera <- SCHEMA_variants[SCHEMA_variants$group %in% c("EUR-N (exomes)"),]
EUR_agg_NonNextera <- aggregate(cbind(EUR_NonNextera$ac_case,EUR_NonNextera$ac_ctrl), by = list(EUR_NonNextera$locus),FUN = sum)
EUR_agg_NonNextera <- EUR_agg_NonNextera[EUR_agg_NonNextera$V1 + EUR_agg_NonNextera$V2 <= 5,]
EUR_agg_NonNextera$gene_id = EUR_NonNextera$gene_id[match(EUR_agg_NonNextera$Group.1,EUR_NonNextera$locus)]
EUR_agg_NonNextera$consequence = EUR_NonNextera$consequence[match(EUR_agg_NonNextera$Group.1,EUR_NonNextera$locus)]
names(EUR_agg_NonNextera) <- c("locus","ac_case","ac_ctrl","gene_id","consequence")

#PTV wrangling
scz_ptv_sumstats_Nextera <- wrangle_sumstats(EUR_agg_Nextera,
                                             8874,
                                             19074+23561,
                                             EUR_agg_Nextera$consequence %in% c("stop_gained",
                                                                                "frameshift_variant",
                                                                                "splice_acceptor_variant",
                                                                                "splice_donor_variant"))
scz_ptv_sumstats_NonNextera <- wrangle_sumstats(EUR_agg_NonNextera,
                                                    7277,
                                                    11187,
                                                    EUR_agg_NonNextera$consequence %in% c("stop_gained",
                                                                                       "frameshift_variant",
                                                                                       "splice_acceptor_variant",
                                                                                       "splice_donor_variant"))

#Missense wrangling
scz_missenseMPC2_sumstats_Nextera <- wrangle_sumstats(EUR_agg_Nextera,
                                             8874,
                                             19074+23561,
                                             EUR_agg_Nextera$consequence %in% c("missense_variant_mpc_2-3",
                                                                                "missense_variant_mpc_>=3"))
scz_missenseMPC2_sumstats_NonNextera <- wrangle_sumstats(EUR_agg_NonNextera,
                                                7277,
                                                11187,
                                                EUR_agg_NonNextera$consequence %in% c("missense_variant_mpc_2-3",
                                                                                      "missense_variant_mpc_>=3"))

#Synonymous wrangling
scz_synonymous_sumstats_Nextera <- wrangle_sumstats(EUR_agg_Nextera,
                                                      8874,
                                                      19074+23561,
                                                      EUR_agg_Nextera$consequence %in% c("synonymous_variant"))
scz_synonymous_sumstats_NonNextera <- wrangle_sumstats(EUR_agg_NonNextera,
                                                         7277,
                                                         11187,
                                                         EUR_agg_NonNextera$consequence %in% c("synonymous_variant"))



#PTV BHR Models
scz_ptv_bhr_Nextera <- bhr::BHR(scz_ptv_sumstats_Nextera, 
                           annotations = list(baseline_model),
                           num_blocks = 100,
                           mode = "univariate", 
                           output_jackknife_h2 = TRUE,
                           all_models = TRUE)

scz_ptv_bhr_NonNextera <- bhr::BHR(scz_ptv_sumstats_NonNextera, 
                                   annotations = list(baseline_model),
                                   num_blocks = 100,
                                   mode = "univariate", 
                                   output_jackknife_h2 = TRUE,
                                   all_models = TRUE)

#Missense BHR Models
scz_missenseMPC2_bhr_Nextera <- bhr::BHR(scz_missenseMPC2_sumstats_Nextera, 
                                annotations = list(baseline_model),
                                num_blocks = 100,
                                mode = "univariate", 
                                output_jackknife_h2 = TRUE,
                                all_models = TRUE)

scz_missenseMPC2_bhr_NonNextera <- bhr::BHR(scz_missenseMPC2_sumstats_NonNextera, 
                                   annotations = list(baseline_model),
                                   num_blocks = 100,
                                   mode = "univariate", 
                                   output_jackknife_h2 = TRUE,
                                   all_models = TRUE)

#Synonymous BHR Models
scz_synonymous_bhr_Nextera <- bhr::BHR(scz_synonymous_sumstats_Nextera, 
                                annotations = list(baseline_model),
                                num_blocks = 100,
                                mode = "univariate", 
                                output_jackknife_h2 = TRUE,
                                all_models = TRUE)

scz_synonymous_bhr_NonNextera <- bhr::BHR(scz_synonymous_sumstats_NonNextera, 
                                   annotations = list(baseline_model),
                                   num_blocks = 100,
                                   mode = "univariate", 
                                   output_jackknife_h2 = TRUE,
                                   all_models = TRUE)

#Functions for meta-analysis
get_metabeta <- function(betas, ses){
  we = 1 / (ses)^2
  return(sum(betas * we) / sum(we))
}
get_metase <- function(ses){
  we = 1 / (ses)^2
  return(sqrt(1/sum(we)))
}

#Function for meta-analyzing nextera and non-nextera samples, and extracting key parameters of interest

SCZ_meta_analysis <- function(modelNextera,modelNonNextera,sumstatsNextera,sumstatsNonNextera, liability = TRUE){
  

  #calculate factor for conversion to liability scale.
  #Population prevalence is from Charlson et al, 2016, Schizophrenia Bulletin
  nextera_prevalence = 8874/(19074+23561+8874)
  nonnextera_prevalence = 7277/(11187+7277)
  population_prevalence = 0.0028
  
  obs2lia_factor <- function(K, P){
    X <- qnorm(K,lower.tail=FALSE)
    z <- (1/sqrt(2*pi))*(exp(-(X**2)/2))
    factor <- (K*(1-K)*K*(1-K))/(P*(1-P)*(z**2))
    return(factor)
  }
  
  nextera_scalingfactor = obs2lia_factor(population_prevalence,nextera_prevalence)
  nonnextera_scalingfactor = obs2lia_factor(population_prevalence,nonnextera_prevalence)
  
  #If observed scale desired, set scaling factor to 1
  if (liability == FALSE){
    nextera_scalingfactor = 1
    nonnextera_scalingfactor = 1
    
  }
  
  #extract heritability from constrained genes in nextera and non-nextera samples
  h2_constrained = c(modelNextera$mixed_model$heritabilities[1,1]*nextera_scalingfactor,
                     modelNonNextera$mixed_model$heritabilities[1,1]*nonnextera_scalingfactor)
  h2_constrained_se = c(modelNextera$mixed_model$heritabilities[2,1]*nextera_scalingfactor,
                        modelNonNextera$mixed_model$heritabilities[2,1]*nonnextera_scalingfactor)
  
  #extract overall heritability in nextera and non-nextera samples
  h2_all = c(modelNextera$mixed_model$heritabilities[1,5]*nextera_scalingfactor,
             modelNonNextera$mixed_model$heritabilities[1,5]*nonnextera_scalingfactor)
  h2_all_se = c(modelNextera$mixed_model$heritabilities[2,5]*nextera_scalingfactor,
                modelNonNextera$mixed_model$heritabilities[2,5]*nonnextera_scalingfactor)
  
  #meta-analyze constrained heritability ("numerator" of fraction explained by constrained genes)
  numerator = get_metabeta(h2_constrained,h2_constrained_se)
  numerator_se = get_metase(h2_constrained_se)
  
  #meta-analyze total heritability ("denominator" of fraction explained by constrained genes)
  denominator = get_metabeta(h2_all,h2_all_se)
  denominator_se = get_metase(h2_all_se)
  
  #get point estimate of fraction explained by constrained genes
  fraction = numerator/denominator
  
  #compute fraction of alleles in constrained annotation
  sumstatsNextera_annot = merge(sumstatsNextera, baseline_model, by.x = "gene", by.y = "gene")
  sumstatsNonNextera_annot = merge(sumstatsNonNextera, baseline_model, by.x = "gene", by.y = "gene")
  sumstatsNextera_annot$burden_score = sapply(sumstatsNextera_annot$variant_variances, function(x) sum(x))
  sumstatsNonNextera_annot$burden_score = sapply(sumstatsNonNextera_annot$variant_variances, function(x) sum(x))
  
  total_variants = sum(sumstatsNextera_annot$burden_score) + sum(sumstatsNonNextera_annot$burden_score)
  constrained_variants = sum(sumstatsNextera_annot$burden_score * sumstatsNextera_annot$baseline_oe1_total5) + sum(sumstatsNonNextera_annot$burden_score*sumstatsNonNextera_annot$baseline_oe1_total5)
  fraction_constrained_variants = constrained_variants/total_variants
  
  #compute constraint enrichment point estimate
  enrichment = fraction/fraction_constrained_variants
  
  #get SE of fraction with delta method
  
  #First, compute covariance matrix of (h2 total, h2 constrained), for nextera and non-nextera
  sigma = matrix(data = NA, nrow = 2,ncol = 2)
  sigma[1,1] = numerator_se^2
  sigma[2,2] = denominator_se^2
  
  #jackknife variances of h2 constrained and h2 total, and jackknfie covariance between h2 constrained and h2 total
  Nextera_constrained_jackknife =  modelNextera$subthreshold_genes$jackknife_h2[1,]*nextera_scalingfactor
  Nextera_total_jackknife =  modelNextera$subthreshold_genes$jackknife_h2[5,]*nextera_scalingfactor
  
  variance_constrained_Nextera = ((length(Nextera_constrained_jackknife) -1)/length(Nextera_constrained_jackknife))*sum((Nextera_constrained_jackknife - mean(Nextera_constrained_jackknife))^2)
  variance_total_Nextera = ((length(Nextera_total_jackknife) -1)/length(Nextera_total_jackknife))*sum((Nextera_total_jackknife - mean(Nextera_total_jackknife))^2)
  covariance_Nextera_jackknife = ((length(Nextera_constrained_jackknife) -1)/length(Nextera_constrained_jackknife)) * sum((Nextera_constrained_jackknife - mean(Nextera_constrained_jackknife))*(Nextera_total_jackknife - mean(Nextera_total_jackknife)))
  
  sigma_nextera = matrix(data = NA, nrow = 2,ncol = 2)
  sigma_nextera[1,1] = variance_constrained_Nextera
  sigma_nextera[2,2] = variance_total_Nextera
  sigma_nextera[1,2] = covariance_Nextera_jackknife
  sigma_nextera[2,1] = covariance_Nextera_jackknife
  
  NonNextera_constrained_jackknife =  modelNonNextera$subthreshold_genes$jackknife_h2[1,]*nonnextera_scalingfactor
  NonNextera_total_jackknife =  modelNonNextera$subthreshold_genes$jackknife_h2[5,]*nonnextera_scalingfactor
  variance_constrained_NonNextera = ((length(NonNextera_constrained_jackknife) -1)/length(NonNextera_constrained_jackknife))*sum((NonNextera_constrained_jackknife - mean(NonNextera_constrained_jackknife))^2)
  variance_total_NonNextera = ((length(NonNextera_total_jackknife) -1)/length(NonNextera_total_jackknife))*sum((NonNextera_total_jackknife - mean(NonNextera_total_jackknife))^2)
  covariance_NonNextera_jackknife = ((length(NonNextera_constrained_jackknife) -1)/length(NonNextera_constrained_jackknife)) * sum((NonNextera_constrained_jackknife - mean(NonNextera_constrained_jackknife))*(NonNextera_total_jackknife - mean(NonNextera_total_jackknife)))
  
  sigma_NonNextera = matrix(data = NA, nrow = 2,ncol = 2)
  sigma_NonNextera[1,1] = variance_constrained_NonNextera
  sigma_NonNextera[2,2] = variance_total_NonNextera
  sigma_NonNextera[1,2] = covariance_NonNextera_jackknife
  sigma_NonNextera[2,1] = covariance_NonNextera_jackknife
  
  #meta-analyze covariance matrices
  sigma_meta = pracma::inv(pracma::inv(sigma_nextera) + pracma::inv(sigma_NonNextera))
  
  
  #get gradient
  dg_dh2c = 1/denominator
  dg_dh2all = (-numerator)/(denominator^2)
  
  fraction_gradient = matrix(c(dg_dh2c,dg_dh2all), ncol = 1)
  
  #compute variance of fraction of heritability explained by constrained genes
  fraction_se = sqrt(t(fraction_gradient) %*% sigma_meta %*% fraction_gradient)
  
  #convert fraction to enrichment
  enrichment_se = fraction_se/fraction_constrained_variants
  
  return(list(bhr_h2 = denominator,
              bhr_h2_se = denominator_se,
              fraction_constrained = fraction,
              fraction_se = fraction_se,
              enrichment_constrained = enrichment,
              enrichment_constrained_se = enrichment_se))
  
}

scz_ptv_bhr_output = SCZ_meta_analysis(scz_ptv_bhr_Nextera,
                                            scz_ptv_bhr_NonNextera,
                                            scz_ptv_sumstats_Nextera,
                                            scz_ptv_sumstats_NonNextera)

scz_missense_bhr_output = SCZ_meta_analysis(scz_missenseMPC2_bhr_Nextera,
                                            scz_missenseMPC2_bhr_NonNextera,
                                            scz_missenseMPC2_sumstats_Nextera,
                                            scz_missenseMPC2_sumstats_NonNextera)

scz_synonymous_bhr_output = SCZ_meta_analysis(scz_synonymous_bhr_Nextera,
                                            scz_synonymous_bhr_NonNextera,
                                            scz_synonymous_sumstats_Nextera,
                                            scz_synonymous_sumstats_NonNextera)

#Gather output
output_df_SCZBP <- data.frame(dx =c(rep("SCZ",3),rep("BP",3)),
                              class = c(c("pLoF","Missense (MPC >2)", "Syn"),c("pLoF","Missense (MPC >2)", "Syn")),
                              bhr_h2 = c(scz_ptv_bhr_output$bhr_h2,
                                         scz_missense_bhr_output$bhr_h2,
                                         scz_synonymous_bhr_output$bhr_h2,
                                         bp_ptv_bhr$mixed_model$heritabilities[1,5]*bp_scalingfactor,
                                         bp_missense_bhr$mixed_model$heritabilities[1,5]*bp_scalingfactor,
                                         bp_synonymous_bhr$mixed_model$heritabilities[1,5]*bp_scalingfactor),
                              bhr_h2_se = c(scz_ptv_bhr_output$bhr_h2_se,
                                            scz_missense_bhr_output$bhr_h2_se,
                                            scz_synonymous_bhr_output$bhr_h2_se,
                                            bp_ptv_bhr$mixed_model$heritabilities[2,5]*bp_scalingfactor,
                                            bp_missense_bhr$mixed_model$heritabilities[2,5]*bp_scalingfactor,
                                            bp_synonymous_bhr$mixed_model$heritabilities[2,5]*bp_scalingfactor),
                              fraction_constrained = c(scz_ptv_bhr_output$fraction_constrained,
                                                       scz_missense_bhr_output$fraction_constrained,
                                                       scz_synonymous_bhr_output$fraction_constrained,
                                                       bp_ptv_bhr$mixed_model$fractions[1,1],
                                                       bp_missense_bhr$mixed_model$fractions[1,1],
                                                       bp_synonymous_bhr$mixed_model$fractions[1,1]),
                              fraction_constrained_se = c(scz_ptv_bhr_output$fraction_constrained_se,
                                                          scz_missense_bhr_output$fraction_constrained_se,
                                                          scz_synonymous_bhr_output$fraction_constrained_se,
                                                          bp_ptv_bhr$mixed_model$fractions[2,1],
                                                          bp_missense_bhr$mixed_model$fractions[2,1],
                                                          bp_synonymous_bhr$mixed_model$fractions[2,1]),
                              enrichment_constrained = c(scz_ptv_bhr_output$enrichment_constrained,
                                                         scz_missense_bhr_output$enrichment_constrained,
                                                         scz_synonymous_bhr_output$enrichment_constrained,
                                                         bp_ptv_bhr$mixed_model$enrichments[1,1],
                                                         bp_missense_bhr$mixed_model$enrichments[1,1],
                                                         bp_synonymous_bhr$mixed_model$enrichments[1,1]),
                              enrichment_constrained_se = c(scz_ptv_bhr_output$enrichment_constrained_se,
                                                            scz_missense_bhr_output$enrichment_constrained_se,
                                                            scz_synonymous_bhr_output$enrichment_constrained_se,
                                                            bp_ptv_bhr$mixed_model$enrichments[2,1],
                                                            bp_missense_bhr$mixed_model$enrichments[2,1],
                                                            bp_synonymous_bhr$mixed_model$enrichments[2,1]))

#fractions of SCZ h2 explained by sig genes
#Load in gnomad information to match gene ids to gene names
gnomad_information <- data.frame(data.table::fread("~/Documents/oconnor_rotation/rarevariantproject/final_manuscript_repo/rv_h2/reference_files/gnomad.v2.1.1.lof_metrics.by_gene.txt",
                                       select = c("gene","gene_id", "chromosome", "start_position",	"end_position", "pLI")))
gnomad_information = gnomad_information[gnomad_information$chromosome != "X",]
gnomad_information = gnomad_information[gnomad_information$chromosome != "Y",]
gnomad_information_counts = data.frame(table(gnomad_information$gene))
gnomad_information_counts = gnomad_information_counts[gnomad_information_counts$Freq == 1,]
gene_information <- gnomad_information[gnomad_information$gene %in% gnomad_information_counts$Var1,]
gene_information$midpoint <- (gene_information$start_position + gene_information$end_position)/2

#9 autosomal SCHEMA genes
SCHEMAgenes = c("SETD1A", "CUL1", "XPO7","TRIO","CACNA1G","SP4","GRIA3","GRIN2A","HERC1","RB1CC1")
SCHEMAgenes_id = gene_information$gene_id[match(SCHEMAgenes,gene_information$gene)]
SCHEMAgenes_id = SCHEMAgenes_id[!is.na(SCHEMAgenes_id)]

#Run BHR with SCHEMA significant genes as fixed effects.
nextera_model_SCHEMA = bhr::BHR(scz_ptv_sumstats_Nextera, 
                           annotations = list(baseline_model),
                           fixed_genes =SCHEMAgenes_id, 
                           num_blocks = 100,
                           mode = "univariate")

nonnextera_model_SCHEMA = bhr::BHR(scz_ptv_sumstats_NonNextera, 
                              annotations = list(baseline_model),
                              fixed_genes =SCHEMAgenes_id, 
                              num_blocks = 100,
                              mode = "univariate")

#Meta analyze fraction of heritability explained by sig genes
meta_frac_SCHEMA = get_metabeta(c(nextera_model_SCHEMA$significant_genes$fraction_burdenh2_significant,
                                  nonnextera_model_SCHEMA$significant_genes$fraction_burdenh2_significant),
                                c(nextera_model_SCHEMA$significant_genes$fraction_burdenh2_significant_se,
                                  nonnextera_model_SCHEMA$significant_genes$fraction_burdenh2_significant_se))

meta_frac_SCHEMA_se = get_metase(c(nextera_model_SCHEMA$significant_genes$fraction_burdenh2_significant_se,
                                   nonnextera_model_SCHEMA$significant_genes$fraction_burdenh2_significant_se))

#rg between scz and bp
nextera_bp_rg = bhr::BHR(mode = "bivariate",
                    trait1_sumstats = scz_ptv_sumstats_Nextera,
                    trait2_sumstats = bp_sumstats_ptv,
                    annotations = list(baseline_model),
                    num_blocks = 100)

nonnextera_bp_rg = bhr::BHR(mode = "bivariate",
                       trait1_sumstats = scz_ptv_sumstats_NonNextera,
                       trait2_sumstats = bp_sumstats_ptv,
                       annotations = list(baseline_model),
                       num_blocks = 100)

write.csv(output_df_SCZBP,"~/Documents/SCZ_BP_output.csv", quote = FALSE, row.names = FALSE)



