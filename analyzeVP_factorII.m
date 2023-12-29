% Load the VP for factor II and analyze results
clear all;

% File where VP are saved
fsave = "28-Dec-2023_FactorII_VP_n-10000_dsgbias_note-dsgshuffle.mat";
%fsave_info = "28-Dec-2023_FactorII_VP_n-10000_note-factorII_info.mat";

factor = 'II';

% Load VPs
VP = load(fsave);
%VP_info = load(fsave_info);
noOC = VP.samplesNoOC;
dsg = VP.samplesDsg;
lev = VP.samplesLev;

% compute differences
diff_lev = lev - noOC;
diff_dsg = dsg - noOC;

% check directions of lev and dsg
pos_diff_dsg = diff_dsg>0;
pos_diff_lev = diff_lev>0;

ids = find(pos_diff_dsg ~= pos_diff_lev);


%% 
% figure specs
w_bin2 = 2;
cmap = parula(5);
c_h = 1; c_kdf = 3; c_samp = 5;
ms = 20;
temp = gray(3);
cgray = temp(2,:);
% plot differences 
figure(13)
yrange = [0,0.15];
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

% counts
figure(14)
clf; 
hold on
histogram(diff_lev, ...
                'BinWidth', w_bin2, 'FaceColor', cmap(2,:))
histogram(diff_dsg, ...
                'BinWidth', w_bin2, 'FaceColor', cmap(4,:))

xlabel(strcat('Factor ', factor,' level difference after OC'))
ylabel('count')
legend('Lev', 'Dsg')
hold off

%% 
% plot pairs
figure(6)
clf;
hold on
ax = gca;
set(ax,'FontSize',18)
ylim([50,180])
xlim([50,180])
% x = y line
temp = xlim(gca);
x = linspace(temp(1),temp(2));
plot(x,x,'color',cgray,'linewidth',2)
plot(noOC, lev, 'linestyle', 'none', ...
    'marker', '.', 'markersize', ms, ...
    'color', cmap(2,:))
plot(noOC, dsg, 'linestyle', 'none', ...
    'marker', '.', 'markersize', ms, ...
    'color', cmap(4,:))
xlabel(strcat('Factor ', factor,' before OC'))
ylabel(strcat('Factor ',factor,' after OC'))
title({'VP pairs', ['Factor ', factor]})
legend('','Lev','Dsg','location','southeast') 
hold off