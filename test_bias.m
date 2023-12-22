clear all;

mu = 0; % location parameter
sigma = 0.001 * 1e4; % scale parameter
N = 1000; % number of samples
pd = makedist('Normal', 'mu', mu, 'sigma', sigma);
samples = round(random(pd, N, 1));

figure(1)
clf;
histogram(samples)

