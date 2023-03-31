import hail as hl

def bhr_genebass_variant_to_gene_lof_syn(var_type, upper_bound, lower_bound, name):
    #Load variants and phenotype files
    genebass_variant = hl.read_matrix_table('gs://ukbb-exome-public/500k/results/variant_results.mt')
    genebass_variant = genebass_variant.rename({'AF.Cases': 'AF_Cases', 'AF.Controls': 'AF_Controls'})
    genebass_variant = genebass_variant.key_cols_by(genebass_variant.phenocode, genebass_variant.coding_description)
    bhr_phenotypes = hl.import_table('~/ms_phenotype_sheet.txt', #This sheet is Supplementary Table 4 from the manuscript
                                types = {'n_cases':hl.tint32, 'n_controls':hl.tint32, 'n_eff':hl.tint32})
    bhr_phenotypes = bhr_phenotypes.key_by(bhr_phenotypes.phenocode, bhr_phenotypes.coding_description)
    genebass_variant = genebass_variant.semi_join_cols(bhr_phenotypes)
    genebass_variant = genebass_variant.annotate_cols(**bhr_phenotypes[genebass_variant.phenocode, genebass_variant.coding_description])
    genebass_variant = genebass_variant.filter_cols(genebass_variant.phenotype_key == "50NA", keep = True)
    
    #Filter variants and count
    genebass_variant = genebass_variant.filter_rows((genebass_variant.call_stats.AF < upper_bound) & (genebass_variant.call_stats.AF >= lower_bound), keep=True)
    genebass_variant = genebass_variant.filter_rows(genebass_variant.annotation == var_type)
    nvar = genebass_variant.count()[0]
    
    #Add BHR-specific parameters
    genebass_variant = genebass_variant.annotate_entries(af_overall = ((genebass_variant.n_cases*genebass_variant.AF_Cases) + (genebass_variant.n_controls*genebass_variant.AF_Controls))/(genebass_variant.n_cases + genebass_variant.n_controls),
                    prevalence = genebass_variant.n_cases/(genebass_variant.n_cases + genebass_variant.n_controls))
    genebass_variant = genebass_variant.annotate_entries(beta_binary = ((2*genebass_variant.prevalence*(genebass_variant.AF_Cases - genebass_variant.af_overall))/hl.sqrt(2*genebass_variant.af_overall*(1-genebass_variant.af_overall)*genebass_variant.prevalence*(1-genebass_variant.prevalence))))
    
    genebass_variant = genebass_variant.annotate_entries(variant_variance = hl.if_else(genebass_variant.trait_type == "continuous",
                                                            ((1/(genebass_variant.SE*hl.sqrt(genebass_variant.n_eff)))**2),
                                                            2*genebass_variant.af_overall*(1-genebass_variant.af_overall)))
                                                                                                                                                    
    genebass_variant = genebass_variant.annotate_entries(beta_per_allele = hl.if_else(genebass_variant.trait_type == "continuous",
                                                          genebass_variant.BETA,
                                                          genebass_variant.beta_binary/(hl.sqrt(genebass_variant.variant_variance))))
                                                                                                                                                               
                 
    #Export gene summary statistics file
    genebass_variant.entries().export('~/dec_bhr_ms_variant_ss_400k_final_thin_withnullburden_'+str(var_type)+'_nvar'+str(nvar)+'_low'+str(lower_bound)+'_high'+str(upper_bound)+'_'+str(name)+'.txt.bgz')
     
     
#variant frequency and function filters
upper = [1e-5]
lower = [0]
names = ['group1']
var_type = ['pLoF']

for variant_type in range(len(var_type)):
    for grouping in range(len(names)):
        bhr_genebass_variant_to_gene_lof_syn(var_type[variant_type], upper[grouping], lower[grouping], names[grouping])
