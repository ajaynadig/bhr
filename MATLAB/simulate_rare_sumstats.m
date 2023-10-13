function [genes, variants, X, yy, popn] = ...
    simulate_rare_sumstats(nn,gg,mm_per_gene,...
    sigmasqSupport,varargin)
%simulate_rare_sumstats simulates rare variant summary statistics
%(gene-level and variant-level) with negative selection and no LD, and true
%effect sizes drawn from a mixture of normal distributions
%
%   Output arguments
%   genes: gene-level summary statistics as a structure with fields:
%       burdenStatistic: sample correlation between Y and gene minor allele
%           burden
%       burdenEffect: true value of gamma (i.e., population correlation)
%       burdenScore: burden score (roughly the combined AF)
%       constraint: observed / expected number of variants in gene
%       noVariants: maximum number of variants (denominator of constraint)
%       noVariantsThreshold: number of variants within MAF bin if specified
%       pqbar: mean heterozygosity of variants in gene within MAF bin
%
%   variants: variant-level summary statistics as a structure.
%
%   Required input arguments:
%   nn: association study sample size
%   gg: number of genes
%   mm_per_gene: max number of variants per gene (number of polymorphic
%   sites will be smaller)
%   sigmasqSupport: possible values for sigmasq, the per-gene effect size
%   variance, as a vector. For example, input "0" to produce null
%   simulations and "[0 1 10]" to produce simulations with a mixture of
%   null genes, small-effect genes, and large-effect genes with effect size
%   variants 10x that of the small-effect genes. Specify optional h2Target
%   to avoid having to put these values on an appropriate scale.
%
%   Optional input arguments:
%   h2Target: desired heritability for variants in the specified AF range
%   maxAF: maximum AF of variants to report sumstats for (effect sizes will
%   also be drawn for other variants but won't be used to compute gene
%   burden statistics)
%   minAF: minimum AF of variants to report sumstats for
%   sigmasqPrior: mixture weight for each element of sigmasqSupport. Gene
%   effect size distribution will be
%   sigmasqPrior(1) * N(0,sigmasqSupport(1)) + ...
%   overdispSupport: variance of overdispersion effects for genes drawn
%   from the same mixture as the gene effects. Effect of variants within
%   each gene will be
%   beta|k ~ N(gamma, overdispSupport(k) * eye)
%   where gamma|k ~ N(0,sigmasqPrior(k)) and eye is the identity matrix.
%   annot: gene-level annotation, allowing for different mixture weights in
%   different annotations. Number of columns of annot should be equal to
%   the number of rows of sigmasqPrior. Number of rows should be gg.
%   popStratVar: variance of population stratification effects. For
%   example, popStratVar=0.1/nn will produce a mean burden chi^2 statistic
%   of 1.1 in null simulations
%   popStratMean: genome-wide mean population stratification effect.
%   Stratification effects are added to each variant, in per-s.d. units,
%   with size alpha ~ N(popStratMean, popStratVar * eye)
%   meanNs: selection strength, as quantified by the mean selection
%   coefficient times the effective population size. A variant with Ns =
%   50, for example, will never reach an allele frequency greater than
%   1/4Ns = 0.005 (see nonneutral_af.m)
%   selectionModel: there are two options, 'stabilizing' and 'directional',
%   for the relationship between effect sizes and selection coefficients
%   noTraits: the number of traits under selection. These are given
%   independent genetic architectures, and results are only reported for
%   one of them, but all of them affect the selection coeffi gcients equally;
%   see Simons et al. 2018 PloS Bio
%   nn_AFS: sample size for calculating allele frequency spectrum
%   (recommended to use default value, which is the study sample size)

p = inputParser;
addRequired(p, 'nn', @isscalar)
addRequired(p, 'gg', @isscalar)
addRequired(p, 'mm_per_gene', @(x)isscalar(x) || all(size(x)==[gg,1]))
addRequired(p, 'sigmasqSupport', @isvector)
addOptional(p, 'h2Target', [], @(x)isscalar(x) || isempty(x))
addOptional(p, 'maxAF', 1, @isscalar)
addOptional(p, 'minAF', 0, @isscalar)
addOptional(p, 'sigmasqPrior', 1, @(x)all(sum(x,2)==1) && size(x,2)==length(sigmasqSupport))
addOptional(p, 'overdispSupport', zeros(size(sigmasqSupport)), @(x)all(size(x)==size(sigmasqSupport)))
addOptional(p, 'annot', ones(gg,1), @(x)size(x,1)==gg)
addOptional(p, 'popStratVar', 0, @isscalar)
addOptional(p, 'popStratMean', 0, @isscalar)
addOptional(p, 'meanNs', 0, @isscalar)
addOptional(p, 'selectionModel', 'stabilizing', @isstr)
addOptional(p, 'noTraits', 1, @(x)isscalar(x) && x>=1)
addOptional(p, 'nn_AFS', nn, @isscalar)
addOptional(p, 'migration_graph', [], @ismatrix)
addOptional(p, 'deme_population_size', [], @iscolumn)
addOptional(p, 'no_generations', [], @isscalar)
addOptional(p, 'mutation_rate', [], @isscalar)
addOptional(p, 'deme_mean_phenotype', [], @iscolumn)

parse(p, nn, gg, mm_per_gene, sigmasqSupport, varargin{:});

% turns p.Results.x into just x
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end

noAnnot = size(annot,2);
if size(sigmasqPrior) ~= [noAnnot, length(sigmasqSupport)]
    error('prior for sigmasq should have size number of annotations x number of components')
end

% allele frequencies within each gene
if isscalar(mm_per_gene)
    mm_per_gene = mm_per_gene * ones(gg,1);
end
mm_tot = sum(mm_per_gene);

% which gene each SNP belongs to
variants.gene = repelem((1:gg)',mm_per_gene);

% Simulate effect sizes for noTraits traits. Only one set of effect sizes
% is kept; others are just used to determine selection coefficients. Traits
% are independent except for shared enrichments if annotations specified
annotSigmasq = zeros(gg,1);
genesigmasq = annotSigmasq;
geneoverdispsigmasq = annotSigmasq;
variants.Ns = 0;
for t = 1:noTraits
    for k = 1:noAnnot
        annotSigmasq(:) = randsample(sigmasqSupport,gg,true,sigmasqPrior(k,:));
        genesigmasq = genesigmasq + annot(:,k) .* annotSigmasq;
        annotSigmasq(:) = randsample(overdispSupport,gg,true,sigmasqPrior(k,:));
        geneoverdispsigmasq = geneoverdispsigmasq + annot(:,k) .* annotSigmasq;
    end

    geneMeanEffect = randn(gg,1) .* sqrt(genesigmasq);

    beta = geneMeanEffect(variants.gene) + ...
        randn(mm_tot,1) .* sqrt(geneoverdispsigmasq(variants.gene));
    if strcmp(p.Results.selectionModel,'stabilizing')
        variants.Ns = variants.Ns + beta.^2;
    elseif strcmp(p.Results.selectionModel,'directional')
        variants.Ns = variants.Ns + max(0,beta);
    else
        error('For selectionModel, choose stabilizing or directional')
    end
end
variants.Ns = variants.Ns * p.Results.meanNs / mean(variants.Ns);

% variant effect sizes
variants.effectPerAllele = beta;

% Simulate allele frequencies directly
variants.AF = nonneutral_af(2 * nn_AFS, variants.Ns);
variants.het = (2 * variants.AF .* (1-variants.AF));

% Simulate allele frequencies using forward simulations
if ~isempty(migration_graph)
    migration_graph = migration_graph ./sum(migration_graph,2);
    noDemes = length(deme_population_size);
    ss = variants.Ns' / sum(deme_population_size);
    AF = repmat(variants.AF',noDemes,1);
    for gen = 1:no_generations
        AF = simulateGeneration(AF,deme_population_size,migration_graph,mutation_rate * ones(size(ss)),ss);
    end
    variants.AF = sum(AF .* deme_population_size)' / sum(deme_population_size);
end

if isempty(migration_graph)
    variants.effect = variants.effectPerAllele .* sqrt(variants.het);

    % add pop strat
    variants.effectEstimate = variants.effect + popStratMean...
        + randn(mm_tot,1) * sqrt(popStratVar);
    
    % add sampling noise
    variants.effectEstimate = variants.effectEstimate + randn(mm_tot,1) / sqrt(nn);
else
    nn_per_deme = mnrnd(nn,deme_population_size/sum(deme_population_size));
    counter = 0;
    X = zeros(nn,mm_tot,'int8');
    yy = zeros(nn,1);
    popn = yy;
    for deme = 1:noDemes
        X(counter+1:counter+nn_per_deme(deme),:) = binornd(...
            2*ones(nn_per_deme(deme),mm_tot),...
            AF(deme,:) .* ones(nn_per_deme(deme),1));
        popn(counter+1:counter+nn_per_deme(deme)) = deme;
        counter = counter + nn_per_deme(deme);
    end
    
    % New allele frequencies
    variants.AF = mean(X)/2;
    variants.het = (2 * variants.AF .* (1-variants.AF))';
    
end

incl = variants.AF <= p.Results.maxAF & variants.AF > p.Results.minAF;
if sum(incl) == 0
    error('No variants sampled in allele frequency bin')
elseif sum(incl) < 1e-3 * mm_tot
    warning('Very few variants sampled in allele frequency bin')
end
variants.effectPerAllele(~incl) = 0;

if ~isempty(migration_graph)
    
    % genetic component
    counter = 0;
    for deme = 1:noDemes
        yy(counter+1:counter+nn_per_deme(deme)) = ...
            single(X(counter+1:counter+nn_per_deme(deme),:))...
            * variants.effectPerAllele;
        counter = counter + nn_per_deme(deme);
    end
    % normalize to mean zero/variance h2
    yy = yy - mean(yy);
    if ~isempty(h2Target)
        yy = yy / std(yy) * sqrt(h2Target);
    end
    
    % add population stratification effects
    counter = 0;
    for deme = 1:noDemes
        yy(counter+1:counter+nn_per_deme(deme)) = ...
            yy(counter+1:counter+nn_per_deme(deme)) + deme_mean_phenotype(deme);
        counter = counter + nn_per_deme(deme);
    end
    yy = yy - mean(yy);
    
    if var(yy) < 1
        yy = yy + randn(nn,1) * sqrt(1 - var(yy));
    else
        error('Y should have variance < 1; check that deme_mean_phenotype makes sense');
    end
    variants.effectEstimate = corr(single(X),yy);
    variants.effectEstimate(isnan(variants.effectEstimate)) = 0;
end

% per-normalized-genotype effects
variants.effect = variants.effectPerAllele .* sqrt(variants.het);

% true h2
h2Burden = sum(geneMeanEffect(variants.gene(incl)).^2 .* variants.het(incl));
if ~isempty(h2Target)
    variants.effect = variants.effect * sqrt(h2Target/h2Burden);
    variants.effectPerAllele = variants.effectPerAllele * sqrt(h2Target/h2Burden);
end

% estimated constraint level for each gene
temp = [0;cumsum(variants.het)];
burdenScoreNothreshold = zeros(gg+1,1);
burdenScoreNothreshold(2:end) = temp(1+cumsum(mm_per_gene));
genes.constraint = (burdenScoreNothreshold(2:end) - burdenScoreNothreshold(1:end-1))...
    ./ mm_per_gene;

% burden scores after thresholding
mm_per_gene_threshold = histcounts(categorical(variants.gene(incl),1:gg));
temp = [0;cumsum(variants.het(incl))];
burdenScore = zeros(gg+1,1);
burdenScore(2:end) = temp(1+cumsum(mm_per_gene_threshold));
genes.burdenScore = burdenScore(2:end) - burdenScore(1:end-1);

% burden statistics
temp = [0;cumsum(sqrt(variants.het(incl)) .*  variants.effectEstimate(incl))];
burdenStatistic = zeros(gg+1,1);
burdenStatistic(2:end) = temp(1+cumsum(mm_per_gene_threshold));
burdenStatistic = burdenStatistic(2:end) - burdenStatistic(1:end-1);
genes.burdenStatistic = burdenStatistic ./sqrt(genes.burdenScore);
genes.burdenStatistic(isnan(genes.burdenStatistic)) = 0;

% true effect sizes
temp = [0;cumsum(sqrt(variants.het(incl)) .*  variants.effect(incl))];
burdenEffect = zeros(gg+1,1);
burdenEffect(2:end) = temp(1+cumsum(mm_per_gene_threshold));
burdenEffect = burdenEffect(2:end) - burdenEffect(1:end-1);
genes.burdenEffect = burdenEffect ./sqrt(genes.burdenScore);
genes.burdenEffect(isnan(genes.burdenEffect)) = 0;

% overdispersion scores
temp = [0;cumsum(variants.het(incl).^2)];
pqbar = zeros(gg+1,1);
pqbar(2:end) = temp(1+cumsum(mm_per_gene_threshold));
pqbar = pqbar(2:end) - pqbar(1:end-1);
genes.pqbar = pqbar ./ genes.burdenScore;
genes.pqbar(isnan(genes.pqbar)) = 0;

genes.noVariants = mm_per_gene;
genes.noVariantsThreshold = mm_per_gene_threshold';


end

