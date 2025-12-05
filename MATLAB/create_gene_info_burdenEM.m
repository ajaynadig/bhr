function gene_info = create_gene_info_burdenEM(num_genes)
%create_gene_info_burdenEM generates a table with dummy gene info for
%burdenEM

gene_info.gene = (1:num_genes)';
gene_info.gene_id = (1:num_genes)';
gene_info.burden_score_uncorrected = ones(num_genes,1);
gene_info.burden_score = ones(num_genes,1);
gene_info.oe_ranking = (1:num_genes)';
gene_info = struct2table(gene_info);
end