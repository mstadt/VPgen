% Testing Factor VIII

clear all;

N = 1e4; % number of samples

mu = 103;
sigma = 4;
pd1 = makedist('Normal', 'mu', mu, 'sigma', sigma);
samples1 = random(pd1,N,1);
mu = 108;
sigma = 23;
pd2 = makedist('Normal', 'mu', mu, 'sigma', sigma);
samples2 = random(pd2,N,1);

figure(1)
clf;
hold on
histogram(samples1)
histogram(samples2)
legend('samples 1', 'samples 2')
hold off
