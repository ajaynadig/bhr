function result = format_genes_burdenEM(genes, prevalence)
%format_genes_burdenEM formats the genes table for burdenEM simulations


if prevalence > 0
    pen = exp(genes.burdenEffectPerAllele) * prevalence;
    pen(genes.burdenScore==0) = 0;
    result.gene_h2 = ((pen-prevalence).^2 .* genes.burdenScore) / ...
                (prevalence * (1-prevalence));
else
    result.gene_h2 = genes.burdenScore .* genes.burdenEffectPerAllele.^2;
end
[~, ranking] = sort(genes.constraint);
result.gene = ranking;
result.effectSize = genes.burdenEffectPerAllele;
result = struct2table(result);
end