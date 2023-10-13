
mm = 100;
demes = 10;
nn = 200 * ones(demes,1);
noGens = 500;

migration = [0.001 0.01 0.1];
noSims = length(migration);

mu = 1e-3 * ones(1,mm);
% ss = exp(-10*rand(1,mm)); ss = sort(ss);
ss = zeros(size(mu));

for sim = 1:noSims
    G = spdiags(ones(demes,3) .* ...
        [migration(sim), 1-2*migration(sim), migration(sim)],-1:1,demes,demes);
    G = G./sum(G,2);
    af = rand(demes,mm);
    for gen = 1:noGens
        af = simulateGeneration(af,nn,G,mu,ss);
    end
    PC = ((1:demes) - (demes+1)/2) * af;
    [~,order] = sort(PC);
    subplot(1,noSims,sim)
    imagesc(af(:,order))
    title(sprintf('migration: %.0e',migration(sim)))
    clim([0 1])
end