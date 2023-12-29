% clear
clear all
%close all

% Generate VP for factor II with biased DSG
% Load NoOC and Lev populations as generated from VP_factorII.m


% set factor
factor = 'II'
note = 'factorII'

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

clearvars -except F_noOC F_lev F_dsg factor note;

%Desired Mean Difference After Treatment (Lev - NoOC):
% Factor II
MEAN_lev = 12; %Table 1 taken from Middeldorp et al. 2000
STD_lev  = 8;  %The Standard Deviation.

MEAN_dsg = 16;
STD_dsg  = 6;
diff_lims_dsg = [MEAN_dsg - 3*STD_dsg, MEAN_dsg + 3*STD_dsg];

N_vp = 1e4; %100 %1e4 %1e4; % how many virtual patients

% shuffle hyperparameters
sigma1 = 0.1 * N_vp; %0.05 * N_vp; %
sigma2 = 0.25 * sigma1; %0.01*sigma1; % variance for high and low values
p = [25, 75]; % percentiles to change bias

MEAN_err = 0.1; % percentage error from given mean
STD_err  = 0.1; % percentage error from given STD



% How many trials do we want at max
MAX_TRIALS = 5e3;

% set random seed
rng(35)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create plots of original distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel density functions
% No OC Factors 
[PDFF_noOC,XF_noOC] = ksdensity(F_noOC);

%Lev Factors
[PDFF_lev,XF_lev] = ksdensity(F_lev);

% Dsg Factors
[PDFF_dsg,XF_dsg] = ksdensity(F_dsg);


figure(1)
clf;
lw = 3;
cmap = parula(5);
c_h = 1; c_kdf = 3; c_samp = 5;
w_bin = 5;
xrange = [40, 190];
yrange = [0,0.075];
subplot(1,3,1)
histogram(F_noOC,'Normalization','pdf', ...
            'BinWidth', w_bin, 'FaceColor', cmap(c_h,:))
hold on
plot(XF_noOC,PDFF_noOC,'color',cmap(c_kdf,:),'linewidth',lw)
xlim(xrange)
ylim(yrange)
temp = strcat('Factor ',factor, ' (No OC)');
title({'Histogram and Kernel Density Function',temp})

subplot(1,3,2)
histogram(F_lev,'Normalization','pdf',...
            'BinWidth', w_bin, 'FaceColor', cmap(c_h,:))
hold on
plot(XF_lev,PDFF_lev,'color',cmap(c_kdf,:),'linewidth',lw)
xlim(xrange)
ylim(yrange)
temp = strcat('Factor ', factor,' (Lev)');
title({'Histogram and Kernel Density Function',temp})

subplot(1,3,3)
histogram(F_dsg,'Normalization','pdf', ...
            'BinWidth', w_bin, 'FaceColor', cmap(c_h,:))
hold on
plot(XF_dsg,PDFF_dsg,'color',cmap(c_kdf,:),'linewidth',lw)
xlim(xrange)
ylim(yrange)
temp = strcat('Factor ',factor,' (Dsg)');
title({'Histogram and Kernel Density Function',temp})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Inverse Kernel Density Functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Range 
pi = linspace(.0001,.9999,10000);

%No OC Factor: 
xiNoOC = ksdensity(F_noOC,pi,'Function','icdf');

%Lev Factors
xiLev = ksdensity(F_lev,pi,'Function','icdf');

% Dsg Factors
xiDsg = ksdensity(F_dsg, pi, 'Function', 'icdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No OC and Lev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsave = "VP/28-Dec-2023_FactorII_VP_n-10000_limitRange_note-factorII.mat"; %"28-Dec-2023_FactorII_VP_n-10000_note-factorII.mat";
VP = load(fsave);
samplesLev = VP.samplesLev;
samplesNoOC = VP.samplesNoOC;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSG with bias to lev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OBJ = 0; % Do we match objective?
NUM_TRIALS = 0; % number of trials
% Pick dsg samples 
xqDsg = rand(1, N_vp);
samplesDsg = sort(interp1(pi,xiDsg,xqDsg)); % ordered dsg samples
while sum(isnan(samplesDsg)) > 0
    % resample if nan value
    %fprintf('nan value, resampling \n')
    xqDsg = rand(1, N_vp);
    samplesDsg = sort(interp1(pi,xiDsg,xqDsg));
end

% compute differences
diff_dsg = samplesDsg - samplesNoOC;

mean_diff_dsg = mean(diff_dsg);
std_dif_dsg = std(diff_dsg);

% Plot starting mean diff and std diff
figure(10);
clf;
subplot(1,2,1)
hold on
c_p = 1;
yline(MEAN_dsg*(1-MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(MEAN_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(MEAN_dsg*(1+MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('mean difference')
title('MEAN')

subplot(1,2,2)
hold on
yline(STD_dsg*(1-STD_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(STD_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(STD_dsg*(1+STD_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('std difference')
title('STD')

% Biased shuffle and redraw
OBJ_MEAN = 0;
last_samp = 0;
diff_lev = samplesLev - samplesNoOC;
% TO DO: bias DSG samples based on what the lev change was for each patient
%diff_lev = samplesLev - samplesNoOC; % difference on lev for each patient
while OBJ~=1 %and(OBJ ~=1, NUM_TRIALS < MAX_TRIALS)
    NUM_TRIALS = NUM_TRIALS + 1;
    if NUM_TRIALS >= MAX_TRIALS
        fprintf('MAX_TRIALS reached \n')
        break;
    end
    % prevent just shuffling if mean is good
    samp_check = NUM_TRIALS - last_samp;

    % if mean is not within tolerance redraw
    if or(OBJ_MEAN == 0, samp_check > 200)
        last_samp = NUM_TRIALS;

        % pick new dsg values (ordered)
        xqDsg = rand(1, N_vp);
        samplesDsg = sort(interp1(pi,xiDsg,xqDsg)); % ordered dsg samples
        while sum(isnan(samplesDsg)) > 0
            % resample if nan value
            %fprintf('nan value, resampling \n')
            xqDsg = rand(1, N_vp);
            samplesDsg = sort(interp1(pi,xiDsg,xqDsg));
        end
    end

    % shuffle the dsg samples (biased)
    %samplesDsg_new = biasedShuffle(samplesDsg, sigma1, sigma2, p);
    samplesDsg_new = biasedShuffle_DSG(samplesDsg, sigma1, sigma2, p, ...
                                        diff_lev, samplesNoOC);
    samplesDsg_new = check_diffs(samplesDsg_new, samplesNoOC,...
                        diff_lims_dsg);

    % compute differences
    diff_dsg_new = samplesDsg_new - samplesNoOC;

    mean_diff = mean(diff_dsg_new);
    std_diff = std(diff_dsg_new);

    % Plot MEAN and STD difference
    figure(10);
    subplot(1,2,1)
    plot(NUM_TRIALS, mean_diff, 'linestyle','none', 'marker', 'o', ...
                                'markersize', 10,...
                                'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    
    subplot(1,2,2)
    plot(NUM_TRIALS, std_diff, 'linestyle','none', 'marker', 'o', 'markersize', 10,...
                                    'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    if mod(NUM_TRIALS,200) == 0
        pause(0.1)
    end

    % Check objectives
    OBJ_MEAN = 0; OBJ_SD = 0;
    if abs((mean_diff - MEAN_dsg)/MEAN_dsg) < MEAN_err
        OBJ_MEAN = 1;
    end
    if abs((std_diff - STD_dsg)/STD_dsg) < STD_err
        OBJ_SD = 1;
    end

    OBJ_test = OBJ_MEAN + OBJ_SD;
    if OBJ_test == 2
        OBJ = 1;
        fprintf('objective reached in %i trials \n', NUM_TRIALS)
    end

    samplesDsg = samplesDsg_new;
end
fprintf('dsg samples done. \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot samples on kernel density function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(1,3,1)
w_bin2 = 2;
histogram(samplesNoOC,'Normalization','pdf', ...
                'BinWidth', w_bin2, 'FaceColor', cmap(c_samp,:))

subplot(1,3,2)
histogram(samplesLev, 'Normalization', 'pdf',...
                'BinWidth', w_bin2, 'FaceColor', cmap(c_samp,:))

subplot(1,3,3)
histogram(samplesDsg, 'Normalization', 'pdf', ...
                'BinWidth', w_bin2, 'FaceColor', cmap(c_samp,:))

legend({'Data', 'Kernel Density function', 'VP'})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot differences for pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(3);

diff_lev = samplesLev - samplesNoOC;
diff_dsg = samplesDsg - samplesNoOC;

yrange = [0,0.15];
clf; 
subplot(1,2,1)
histogram(diff_lev, ...
                'BinWidth', w_bin2, 'FaceColor', cmap(2,:), ...
                'Normalization', 'pdf')

xlabel(strcat('factor ', factor,' level difference after lev'))
ylabel('frequency')
ylim(yrange)
temp = sprintf('LEV \n MEAN DIFF: %0.3f \n STD DIFF: %0.3f',...
        mean(diff_lev), std(diff_lev));
title(temp)

subplot(1,2,2)
histogram(diff_dsg, ...
                'BinWidth', w_bin2, 'FaceColor', cmap(2,:),...
                'Normalization', 'pdf')
xlabel(strcat('Factor ', factor,' level difference after dsg'))
ylabel('frequency')
ylim(yrange)
temp = sprintf('DSG \n MEAN DIFF: %0.3f \n STD DIFF: %0.3f',...
        mean(diff_dsg), std(diff_dsg));
title(temp)

figure(13)
clf; 
hold on
histogram(diff_lev, ...
                'BinWidth', w_bin2, 'FaceColor', cmap(2,:), ...
                'Normalization', 'pdf')
histogram(diff_dsg, ...
                'BinWidth', w_bin2, 'FaceColor', cmap(4,:),...
                'Normalization', 'pdf')

xlabel(strcat('Factor ', factor,' level difference after OC'))
ylabel('density')
legend('Lev', 'Dsg')
ylim(yrange)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(5)
clf;
subplot(1,2,1)
plot(samplesNoOC, samplesLev, 'linestyle', 'none', 'marker', '.', 'markersize', 15)
xlabel(strcat('Factor ', factor,' before OC'))
ylabel(strcat('Factor ',factor,' after Lev'))
title('Lev VP pairs')
ylim([50,180])
xlim([50,180])

subplot(1,2,2)
plot(samplesNoOC, samplesDsg, 'linestyle', 'none', 'marker', '.', 'markersize', 15)
xlabel(strcat('Factor ',factor,' before OC'))
ylabel(strcat('Factor ',factor,' after Dsg'))
title('Dsg VP pairs')
ylim([50,180])
xlim([50,180])

%% pairs on same axes
% plot pairs
figure(6)
clf;
temp = gray(3);
cgray = temp(2,:);
ms = 20;
hold on
ax = gca;
set(ax,'FontSize',18)
ylim([50,180])
xlim([50,180])
% x = y line
temp = xlim(gca);
x = linspace(temp(1),temp(2));
plot(x,x,'color',cgray,'linewidth',2)
plot(samplesNoOC, samplesLev, 'linestyle', 'none', ...
    'marker', '.', 'markersize', ms, ...
    'color', cmap(2,:))
plot(samplesNoOC, samplesDsg, 'linestyle', 'none', ...
    'marker', '.', 'markersize', ms, ...
    'color', cmap(4,:))
xlabel(strcat('Factor ', factor,' before OC'))
ylabel(strcat('Factor ',factor,' after OC'))
title({'VP pairs', ['Factor ', factor]})
legend('','Lev','Dsg','location','southeast') 
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
save_file = 1;
fname = strcat(date, '_Factor', factor, '_VP', '_n-', ...
                num2str(N_vp), '_dsgbias','_note-', note, '.mat');
if isfile(fname)
    save_file = input('file exists. save file? (0/1)');
    if save_file
        note = input('change note. note: ');
        fname = strcat(date, '_Factor', factor, '_VP', '_n-', ...
                    num2str(N_vp), '_dsgbias','_note-', note, '.mat');
    end
end
if save_file
    save(fname, 'samplesNoOC', 'samplesDsg', 'samplesLev')
    fnameinf = strcat(date, '_Factor', factor, '_VP', ...
        '_n-', num2str(N_vp), '_dsgbias', '_note-', note ,'_info.mat');
    save(fnameinf)
    
    fprintf('VP saved to: \n %s \n', fname)
end






