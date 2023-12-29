clear all;

MEAN_lev = 12; %Table 1 taken from Middeldorp et al. 2000
STD_lev  = 8;  %The Standard Deviation.

diff_lims_lev = [MEAN_lev - 3*STD_lev, MEAN_lev + 3*STD_lev];

MEAN_dsg = 16;
STD_dsg  = 6;

diff_lims_dsg = [MEAN_dsg - 3*STD_dsg, MEAN_dsg + 3*STD_dsg];

rng(10)
fsave = "28-Dec-2023_FactorII_VP_n-10000_note-check.mat";
% Load VPs
VP = load(fsave);
%VP_info = load(fsave_info);
noOC = VP.samplesNoOC;
dsg = VP.samplesDsg;
lev = VP.samplesLev;


% Test check_diffs
fprintf('lev \n')
newlev = check_diffs(lev, noOC, diff_lims_lev);
fprintf('dsg \n')
newdsg = check_diffs(dsg, noOC, diff_lims_dsg);