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
    
    % bino(n,p) ~ pois(np)
    i1 = af < bino_threshold;
    [ii,~] = find(i1);
    af(i1) = poissrnd(nn(ii) .* af(i1)) ./ nn(ii);
    
    % bino(n,p) ~ n - pois(n(p-1))
    i2 = af > 1-bino_threshold;
    [ii,~] = find(i2);
    af(i2) = 1 - poissrnd(nn(ii) .* (1-af(i2))) ./ nn(ii);
    
    % bino(n,p) ~ normal(np,np(1-p))
    i3 = min(af,1-af) > bino_threshold;
    [ii,~] = find(i3);
    af(i3) = randn(1,sum(i3(:))) .* ...
        sqrt(af(i3) .* (1-af(i3)) ./ nn(ii))...
        + af(i3);
    af(i3) = max(0,min(1,af(i3)));
end
end