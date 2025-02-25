function [af] = simulateGeneration(af,nn,G,mu,ss)
%simulateGeneration simulates one generation of migration, mutation,
%selection and drift in a set of populations.
% Input arguments:
%   af: allele frequency of each variant in each deme; rows = demes,
%   columns = variants
%   nn: population size of each deme as a column vector
%   G: migration graph, as a demes-by-demes matrix/sparse matrix. Row
%   sums of G should equal 1.
%   mu: mutation rate (forwards = backwards) of each variant as a row
%   vector
%   ss: selection coefficient of each variant as a row vector.

% allele frequency threshold where an approximation to binomial sampling
% will be used instead of slow binomial sampling
bino_threshold = 0.1;

% migration
af = G*af;

% mutation
af = af .* (1-mu) + (1-af) .* mu;

% selection
af = af .* (1-ss) ./ ((1-af) + af .* (1-ss));

% drift
if all(nn < 100)
    af = binornd(nn.*ones(size(mu)),af) ./ nn;
else % approximate binomial sampling with poisson/gaussian
    af = approx_binornd(nn, af, bino_threshold) ./ nn;
end
end