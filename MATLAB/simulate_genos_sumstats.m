% This script simulates individual-level genotypes + phenotypes for the manuscript
% Weiner*, Nadig* et al medRxiv

% clear
cd /broad/oconnor/h2burden/matlab/

% Which replicate number (used for naming sumstats files)
rep = str2num(getenv('SGE_TASK_ID'));
rng(rep)

% This script will produce one gene sumstats file for each simulation,
% where the number of simulations is
% length(h2BurdenTrue) * length(max_af) * length(titles)
save_path = '/broad/oconnor/h2burden/newsims_bugfixed/';

gg = 1e3; % no. genes

% number of traits under pleiotropic selection
noTraitsUnderSelection = 1;

% no. snps per gene
mm_per_gene = randi(1e2,gg,1);

% true h2 burden; run null and non-null simulations
h2BurdenTrue = [0];

% AF bins; run one simulation for each MAF bin
max_af = 1e-3; % maximum MAF
min_af = 0; % minimum AF

% Description of each simulation. Parameters below correspond to each
% simulation.
titles = {'no_popstrat','north_south','north_bottlenecked','local_bottlenecked',...
    'no_popstrat_selection','north_south_selection',...
    'north_bottlenecked_selection','local_bottlenecked_selection'};

% sample size
nn = 2e4 * ones(1,length(titles));

% mean strength of selection (Ns)
selectionStrength = [zeros(1,4) 10 * ones(1,4)];

% Median of gene effect-size distribution
sparsity = 1e1 * ones(1,length(titles));

% Width of gene effect-size distribution (2 narrow, 10 wide)
width = 2 * ones(1,length(titles));

% overdispersion effects
overdisp = zeros(1,length(titles));

% demographic model parameters
noGens = 1000;
mu = 1e-5;

% three-major-deme model with one isolated deme
nnPerDeme_cells = repmat({[1e4 1e4 1e4]',[1e4 1e4 1e4]',[5e3 1e4 2e4]',[1e4 1e4 1e4 1e3]'},1,2);
G_cells = repmat({sparse([.99 .01 0; .01 .98 .01; 0 .01 .99]);
    sparse([.99 .01 0; .01 .98 .01; 0 .01 .99]);
    sparse([.99 .01 0; .01 .98 .01; 0 .01 .99]);
    sparse([.99 .01 0 0; .01 .97 .01 .01; 0 .01 .99 0; 0 .01 0 .99])},2,1);
yy_cells = repmat({[0 0 0]',[0.5 0 -0.5]',[0.5 0 -0.5]',[0 0 0 1]'},1,2);

for nval=1:length(nn) % loop over simulation setups (titles)

    % mixture component variances
    sigmasq = [1/width(nval) 1 width(nval) 0];

    % probability each mixture cpt (cols) in each annotation (rows), with sum
    % at most 1
    prior = 1/sparsity(nval) ./ sigmasq(1:end-1);
    if sum(prior)>1; error('prior should sum to <=1'); end
    prior = [prior, 1 - sum(prior)]; % null component

    for sim = 1:length(h2BurdenTrue) % loop over true h2 values
        filename = [save_path,titles{nval},'.h2_',num2str(h2BurdenTrue(sim)),'.rep_',num2str(rep)];

        if ~isfile([filename,'.genos.csv'])
            % simulate gene/variant level sumstats
            [genes, variants, X, y, popn] = ...
                simulate_rare_sumstats(nn(nval),gg,mm_per_gene,sigmasq,...
                'overdispSupport',overdisp(nval)*sigmasq,...
                'sigmasqPrior',prior,...
                'maxAF',max_af,...
                'minAF',min_af,...
                'h2Target',h2BurdenTrue(sim),...
                'meanNs', selectionStrength(sim),...
                'noTraits',noTraitsUnderSelection,...
                'selectionModel','stabilizing',...
                'migration_graph',G_cells{nval},...
                'deme_population_size',nnPerDeme_cells{nval},...
                'mutation_rate',mu,...
                'deme_mean_phenotype',yy_cells{nval},...
                'no_generations',noGens);

            h2BurdenEst = run_BHR(genes,0:.25:1);
            disp(h2BurdenEst)
            h2BurdenEst_noBins = run_BHR(genes,0:1);
            disp(h2BurdenEst_noBins)

            % save to file
            writetable(struct2table(genes),[filename,'.sumstats.csv']);
            writePhenFile(y,[filename,'.phen']);
            writePhenFile(popn,[filename,'.covar']);
            writematrix(X,[filename,'.genos.csv']);
        end
    end


end

