function rtcurves_stats(sourcedir)
%RTCURVES_STATS compares cue response of cholinergic neurons partitioned to
%different reaction time quartiles.
%
%   RTCURVES_STATS(SOURCEDIR) compares cue response of cholinergic neurons 
%   partitioned to different reaction time quartiles. It loads psth data of 
%   cholinergic neurons from SOURCEDIR and calculates cue response corresponding
%   to the 1st, 2nd, 3rd and 4th quartile of the reaction time districution.
%   Cue responses are compared with a one-way ANOVA test.
%
% See also ANOVA1 and ULTIMATE_PSTH.

%   Panna Hegedus, Balazs Hangya
%   Insititute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   02-09-2021

% Load data
load(fullfile(sourcedir, 'allpsth_data.mat')); % load previously saved data

inx = find(time>0 & time<=0.5); % time indices

% Extract firing rate from PSTH
allspsth_c500 = allspsth(:,:,inx);
allspsth_c500_q1 = mean(allspsth(1,:,inx),3);
allspsth_c500_q2 = mean(allspsth(2,:,inx),3);
allspsth_c500_q3 = mean(allspsth(3,:,inx),3);
allspsth_c500_q4 = mean(allspsth(4,:,inx),3);

mean_cueresp = [mean(allspsth_c500_q1) mean(allspsth_c500_q2) mean(allspsth_c500_q3) mean(allspsth_c500_q4)];
sem_cueresp = [std(allspsth_c500_q1)/5 std(allspsth_c500_q2)/5 std(allspsth_c500_q3)/5 std(allspsth_c500_q4)/5];

% Bar plot
figure; 
bar([1 2 3 4], mean_cueresp);
hold on
errorbar([1 2 3 4], mean_cueresp, sem_cueresp)

% ANOVA
[p, table] = anova1([allspsth_c500_q1' allspsth_c500_q2' allspsth_c500_q3' allspsth_c500_q4'])