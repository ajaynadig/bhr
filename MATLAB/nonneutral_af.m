function [pp] = nonneutral_af(nn,beta,varargin)
%nonneutral_af samples allele frequencies from an approximate non-neutral
%spectrum. To do so, it samples from a neutral spectrum and sets to zero
%any sample frequencies that are above 4Ns.
%   Output arguments:
%   pp: within-sample allele frequencies.
%   
%   Required input arguments:
%   nn: study sample size
%   beta: effect size of each variant, as a variants-by-traits matrix.
%   Selection coefficient of each variant will be determined by taking sum of
%   the squares along the rows of this matrix
% 
%   Optional input arguments:
%   meanNs: mean selection coefficient times effective population size
%   maxAF: variants above this maximum are always discarded

p = inputParser;
addRequired(p, 'nn', @isscalar)
addRequired(p, 'beta', @ismatrix)
addOptional(p, 'meanNs', [], @isscalar)
addOptional(p, 'maxAF', 1,  @isscalar)

parse(p, nn, beta, varargin{:});
maxAF = p.Results.maxAF;
mm = length(beta);

Ns = sum(beta.^2,2);
if ~isempty(p.Results.meanNs)
    Ns = Ns * abs(p.Results.meanNs) / median(Ns);
end
max_af = min(p.Results.maxAF, 1./(4*Ns));

nmax=nn*p.Results.maxAF;
weight = 1./(1:nmax)';
weight = weight/sum(weight);
pp = randsample((1:nmax)'/nn, mm, true, weight);
pp(pp > max_af) = 0;


end

