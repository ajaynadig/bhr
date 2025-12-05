function formatted = format_variants_burdenEM(...
    variants_table, num_samples, num_cases, gene_constraint)
%format_variants_burdenEM formats the simulated variants_table for burdenEM

num_variants = length(variants_table.gene);

% Genes are re-indexed to be sorted by their constraint; this way, a shared
% 'annotation' file can be used across simulation runs where the constraint
% 'annotation' is just 1,2,...,num_genes
[~, gene_ranking] = sort(gene_constraint);
formatted.gene = gene_ranking(variants_table.gene);

formatted.POS = (1:num_variants)';
formatted.CHR = ones(num_variants, 1);
formatted.beta = variants_table.effectEstimate ./ sqrt(variants_table.het);
formatted.N = num_samples * ones(num_variants, 1);
formatted.variant_variance = variants_table.het;
formatted.AF = variants_table.AF;
formatted.trait_type = cell(num_variants, 1);
if num_cases > 0
    formatted.trait_type(:) = {'binary'};
    formatted.AC_cases = 2 * num_cases * variants_table.AF_cases;
    formatted.prevalence = num_cases / num_samples * ones(num_variants, 1);
else
    formatted.trait_type(:) = {'continuous'};
end

is_polymorphic = formatted.variant_variance > 0;
formatted = struct2table(formatted);
formatted = formatted(is_polymorphic,:);

end