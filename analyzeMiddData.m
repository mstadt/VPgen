clear all;

factor = 'VII'

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

F_noOC = [F_noOC1;F_noOC2];

% Factor VII
MEAN_beforeLev = 99;
SD_beforeLev = 40;
MEAN_beforeDsg = 100;
SD_beforeDsg = 17;



MEAN_lev = 12; %Table 1 taken from Middeldorp et al. 2000
STD_lev  = 15;  %The Standard Deviation.

MEAN_dsg = 32;
STD_dsg  = 10;