% This script simulates rare variant summary statistics for the manuscript
% Weiner*, Nadig* et al medRxiv

clear

% Which replicate number (used for naming sumstats files)
rep = 1;

% This script will produce one gene sumstats file for each simulation,
% where the number of simulations is 
% length(h2BurdenTrue) * length(max_af) * length(titles)
save_path = '/path/to/save/sumstats/files/';

gg = 1.8e4; % no. genes

% number of traits under pleiotropic selection
noTraitsUnderSelection = 100;

% no. snps per gene
mm_per_gene = randi(1e3,gg,1);

% true h2 burden; run null and non-null simulations
h2BurdenTrue = [0 0.005];

% AF bins; run one simulation for each MAF bin
max_af = [1e-5 1e-4 1e-3 1e-2]; % maximum MAF
min_af = [0 max_af(1:3)]; % minimum AF

% Width of effect-size distribution (2 narrow, 10 wide)
width = 5 * ones(1,10);

% Description of each simulation. Parameters below correspond to each
% simulation.
titles = {'realistic','small_N','large_N','strong_selection','no_selection',...
    'strong_popstrat','no_popstrat','overdispersion','more_polygenic','less_polygenic'};

% sample size
nn = [5e5 1e5 2e6 5e5 5e5 5e5 5e5 5e5 5e5 5e5]; 

% mean strength of selection (Ns)
selectionStrength = 10 * [1 1 1 10 0 1 1 1 1 1];

% Median of gene heritability distribution
sparsity = [2e3 * ones(1,8) 5e2 1e4];

% overdispersion effects
overdisp = [0 0 0 0 0 0 0 1 0 0]';

% population stratification effects
popstratmean = [1 1 1 1 1 10 0 1 1 1] * 1e-5;
popstratvar = [1 1 1 1 1 10 0 1 1 1] * 1e-7;


clear h2est_noAnnot h2est_Annot h2est_AnnotObs
for nval=1:length(nn) % loop over simulation setups (titles)
    
    % mixture component variances
    sigmasq = [1/width(nval) 1 width(nval) 0];
    
    % probability each mixture cpt (cols) in each annotation (rows), with sum
    % at most 1
    prior = 1/sparsity(nval) ./ sigmasq(1:end-1);
    if sum(prior)>1; error('prior should sum to <=1'); end
    prior = [prior, 1 - sum(prior)]; % null component
    
    for sim = 1:length(h2BurdenTrue) % loop over true h2 values
        for bin = 1:length(max_af) % loop over AF bins
            
            % simulate gene/variant level sumstats
            [genes, variants] = ...
                simulate_rare_sumstats(nn(nval),gg,mm_per_gene,sigmasq,...
                'overdispSupport',overdisp(nval)*sigmasq,...
                'sigmasqPrior',prior,...
                'maxAF',max_af(bin),...
                'minAF',min_af(bin),...
                'meanNs', 0,...
                'popStratVar',popstratvar(nval),...
                'popStratMean',popstratmean(nval),...
                'h2Target',h2BurdenTrue(sim),...
                'meanNs', selectionStrength(nval),...
                'noTraits',noTraitsUnderSelection,...
                'selectionModel','stabilizing');
            
            % save to file
            writetable(struct2table(genes),[save_path,titles{nval},'.MAFbin=',num2str(bin),'.h2=',num2str(h2BurdenTrue(sim)),'.rep=',num2str(rep),'.csv']);
            
        end
    end
    
end

