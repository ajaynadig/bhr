function x = approx_binornd(N, p, bino_threshold)
%APPROX_BINORND performs approximate binomial sampling by using either
% a Poisson or a Normal approximation; Poisson is chosen if p (or 1-p) is
% smaller than the threshold. N and p should have the same shape.
    
if nargin < 3
    bino_threshold = 0.1;
end

x = zeros(size(p));

i1 = p < bino_threshold;
x(i1) = poissrnd(N(i1) .* p(i1));

% bino(N,p) ~ N - pois(N(p-1))
i2 = p > 1-bino_threshold;
x(i2) = N(i2) - poissrnd(N(i2) .* (1-p(i2)));

% bino(n,p) ~ normal(np,np(1-p))
i3 = min(p,1-p) >= bino_threshold;
x(i3) = randn(sum(i3(:)),1) .* ...
    sqrt(p(i3) .* (1-p(i3)) .* N(i3))...
    + p(i3).* N(i3);
x(i3) = max(0,min(N(i3),round(x(i3))));

% small N: just use binomial sampling
i4 = N < 100;
x(i4) = binornd(round(N(i4)), p(i4));

end