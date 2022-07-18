import hail as hl

def bhr_genebass_variant_to_gene_mis(missense_type_list, missense_name, upper_bound, lower_bound, name):
    #Load variants and phenotype files
    genebass_variant = hl.read_matrix_table('gs://ukbb-exome-public/500k/results/variant_results.mt')
    genebass_variant = genebass_variant.rename({'AF.Cases': 'AF_Cases', 'AF.Controls': 'AF_Controls'})
    genebass_variant = genebass_variant.key_cols_by(genebass_variant.phenocode, genebass_variant.coding_description)
    bhr_phenotypes = hl.import_table('gs://ms_phenotype_sheet.txt', #This sheet is Supplementary Table 4 from the manuscript
                                types = {'n_cases':hl.tint32, 'n_controls':hl.tint32, 'n_eff':hl.tint32})
    bhr_phenotypes = bhr_phenotypes.key_by(bhr_phenotypes.phenocode, bhr_phenotypes.coding_description)
    genebass_variant = genebass_variant.semi_join_cols(bhr_phenotypes)
    genebass_variant = genebass_variant.annotate_cols(**bhr_phenotypes[genebass_variant.phenocode, genebass_variant.coding_description])
    genebass_variant = genebass_variant.filter_rows((genebass_variant.call_stats.AF < upper_bound) & (genebass_variant.call_stats.AF >= lower_bound), keep=True)
    
    #Filter variants
    vep_ht = hl.read_table("gs://ukbb-exome-public/500k/results/vep.ht/")
    genebass_variant = genebass_variant.filter_rows(genebass_variant.annotation == "missense")
    genebass_variant = genebass_variant.annotate_rows(vep = vep_ht[genebass_variant.row_key].vep)
    genebass_variant = genebass_variant.filter_rows(hl.literal(missense_type_list).contains(genebass_variant.vep.transcript_consequences.polyphen_prediction[0]))
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
                                                                                                                                                               
    #Aggregate variants into genes
    genebass_gene = (genebass_variant.group_rows_by(genebass_variant.gene)).aggregate(
            variant_variances = hl.agg.collect(genebass_variant.variant_variance),
            betas = hl.agg.collect(genebass_variant.beta_per_allele))
        
    #Export gene summary statistics file
    genebass_gene.entries().export('gs://bhr_ms_gene_ss_400k_final_withnullburden_missense-'+str(missense_name)+'_nvar'+str(nvar)+'_low'+str(lower_bound)+'_high'+str(upper_bound)+'_'+str(name)+'.txt.bgz')
    
#variant frequency and function filters
upper = [1e-5,     1e-4,     1e-3,      1e-2,    1]
lower = [0,        1e-5,     1e-4,      1e-3,    0.05]
names = ['group1', 'group2', 'group3', 'group4', 'common']
missense_functional = [['probably_damaging', 'possibly_damaging'], ['benign']]
missense_names =       ['notbenign',                               'benign']

for missense_type in range(len(missense_names)):
    for grouping in range(len(names)):
        bhr_genebass_variant_to_gene_mis(missense_functional[missense_type], missense_names[missense_type], upper[grouping], lower[grouping], names[grouping])
