clear all; 

% set factor
factor = 'VII'
note = 'factorVII'

rng(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data and Set Ranges %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Data
temp = strcat('data/Factor', factor, '_noOC1.csv');
csvtable1 =  readtable(temp);
temp = strcat('data/Factor', factor, '_lev.csv');
csvtable2 =  readtable(temp);
temp = strcat('data/Factor', factor, '_noOC2.csv');
csvtable3 =  readtable(temp);
temp = strcat('data/Factor', factor, '_dsg.csv');
csvtable4 =  readtable(temp);
F_noOC1 = csvtable1.Var2; 
F_lev   = csvtable2.Var2;
F_noOC2 = csvtable3.Var2;
F_dsg   = csvtable4.Var2;

F_noOC = [F_noOC1; F_noOC2]; % merge noOC data together
allIndices = 1:length(F_noOC);

MEAN_lev = 12; %Table 1 taken from Middeldorp et al. 2000
STD_lev  = 15;  %The Standard Deviation.

MEAN_dsg = 32;
STD_dsg  = 10;

MEAN_err= 0.5; % difference from actual value

rng(1);

fname = '02-Jan-2024_FactorVII_newNoOC_lev_dsg.mat';

% Load best indices
best_inds = load(fname);