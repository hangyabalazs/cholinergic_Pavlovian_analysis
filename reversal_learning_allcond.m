function reversal_learning_allcond(cbname, resdir, isbehav)
%REVERSAL_LEARNING_ALLCOND(CBNAME, RESDIR, ISBEHAV)   Average lick PETH.
%   REVERSAL_LEARNING_ALLCOND plots average lick PETH (beam break time stamps
%   aligned to an event) for all conditions of a probabilistic pavlovian
%   condiitoning task (likely reward ,unlikely punishment, likely punishment
%   unlikely reward). All sessions from mice are used.
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   07-June-2022


choosecb(cbname)
allsessions = listtag('allsessions');
NumSessions = size(allsessions,1);

if ~isfolder(resdir)
    mkdir(resdir)
end

% Time window
wn = [-3 3];   % in seconds
twin = [0 1];
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector
tinx = [find(time==twin(1)) find(time==twin(2))];

% PETH for all 4 conditions
[Hit_allpsth, Hit_allpsth2, FA_allpsth, FA_allpsth2, Hit_spt, FA_spt, Hit_spt2, FA_spt2] = deal([]);
for iS = 1:NumSessions
    sessionid = allsessions(iS,:);
    [spsth_hit, spt_hit, spsth_hit2, spt_hit2, spsth_fa, spt_fa, spsth_fa2, spt_fa2] = main2(sessionid,wn,dt);
    Hit_allpsth = [Hit_allpsth; spsth_hit];
    Hit_allpsth2 = [Hit_allpsth2; spsth_hit2];
    FA_allpsth = [FA_allpsth; spsth_fa];
    FA_allpsth2 = [FA_allpsth2; spsth_fa2];
    
    Hit_spt = [Hit_spt; mean(sum(spt_hit(:,tinx(1):tinx(2)),2))];
    Hit_spt2 = [Hit_spt2; mean(sum(spt_hit2(:,tinx(1):tinx(2)),2))];
    FA_spt = [FA_spt; mean(sum(spt_fa(:,tinx(1):tinx(2)),2))];
    FA_spt2 = [FA_spt2; mean(sum(spt_fa2(:,tinx(1):tinx(2)),2))];
end

boxstat2(Hit_spt, Hit_spt2, 'T1', 'T2', 0.05, 'paired')
boxstat2(Hit_spt, FA_spt, 'T1', 'T3', 0.05, 'paired')
boxstat2(Hit_spt, FA_spt2, 'T1', 'T4', 0.05, 'paired')
boxstat2(Hit_spt2, FA_spt, 'T2', 'T3', 0.05, 'paired')
boxstat2(Hit_spt2, FA_spt2, 'T2', 'T4', 0.05, 'paired')
boxstat2(FA_spt, FA_spt2, 'T3', 'T4', 0.05, 'paired')

figure;
bar([1 2 3 4], [nanmedian(Hit_spt) nanmedian(Hit_spt2) nanmedian(FA_spt) nanmedian(FA_spt2)])
sem1 = se_of_median(Hit_spt);
sem2 = se_of_median(Hit_spt2);
sem3 = se_of_median(FA_spt);
sem4 = se_of_median(FA_spt2);

hold on
errorbar([1 2 3 4], [nanmedian(Hit_spt) nanmedian(Hit_spt2) nanmedian(FA_spt) nanmedian(FA_spt2)], [sem1 sem2 sem3 sem4], [sem1 sem2 sem3 sem4])
fnm = fullfile(resdir,'Barplot_allcond_lickPSTH.fig');
saveas(gcf,fnm)
set(gcf,'renderer','painters')
fnm = fullfile(resdir,'Barplot_allcond_lickPSTH.eps');
saveas(gcf,fnm)

% Plot & save
G = figure;
green = [51 204 51] / 255;   % colors for plotting
blue = [0 0 1];
red = [216 41 0] / 255;
dark_red = [0.5 0 0];
errorshade(time,nanmean(Hit_allpsth),nanstd(Hit_allpsth)/sqrt(size(Hit_allpsth,1)),...
    'LineColor',green,'ShadeColor',green)
hold on
errorshade(time,nanmean(Hit_allpsth2),nanstd(Hit_allpsth)/sqrt(size(Hit_allpsth,1)),...
    'LineColor',blue,'ShadeColor',blue)
hold on
errorshade(time,nanmean(FA_allpsth),nanstd(FA_allpsth)/sqrt(size(FA_allpsth,1)),...
    'LineColor',red,'ShadeColor',red)
hold on
errorshade(time,nanmean(FA_allpsth2),nanstd(FA_allpsth2)/sqrt(size(FA_allpsth2,1)),...
    'LineColor',dark_red,'ShadeColor',dark_red)
set(G,'renderer', 'painters')
fnm = fullfile(resdir,['allanimals_average_lickPSTH_allcond_TT_.fig']);
saveas(G,fnm)
fnm = fullfile(resdir,['allanimals_average_lickPSTH_allcond_TT_.eps']);
saveas(G,fnm)
fnm = fullfile(resdir,['allanimals_average_lickPSTH_allcond_TT_.jpg']);
saveas(G,fnm)

function [spsth_hit, spt_hit, spsth_hit2, spt_hit2, spsth_fa, spt_fa, spsth_fa2, spt_fa2] = main2(cellid,wn,dt)

% Filter input
filterinput_hit = 'TrialType==1&AllReward==1';
filterinput_hit2 = 'TrialType==1&Punishment==1';
filterinput_fa = 'TrialType==2&AllReward==2';
filterinput_fa2 = 'TrialType==2&Punishment==2';


% Calcualte lick PSTH
[~, spsth_hit, ~, ~, spt_hit] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.05,'event_filter','custom','filterinput',filterinput_hit,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_fa, ~,~,spt_fa] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.05,'event_filter','custom','filterinput',filterinput_fa,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_hit2, ~,~,spt_hit2] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.05,'event_filter','custom','filterinput',filterinput_hit2,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_fa2, ~,~,spt_fa2] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.05,'event_filter','custom','filterinput',filterinput_fa2,...
    'isadaptive',0,'maxtrialno',Inf);
