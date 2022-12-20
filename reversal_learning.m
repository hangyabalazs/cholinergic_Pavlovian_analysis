function reversal_learning(cbname, resdir, isbehav)
%REVERSAL_LEARNING(CBNAME, RESDIR, ISBEHAV)   Average lick PETH.
%   REVERSAL_LEARNING plots average lick PETH (beam break time stamps
%   aligned to an event). All sessions from mice are used. ISBEHAV controls
%   the need of preprocessing behavioral data (1=yes, 0=no).
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   02-Jun-2022

% Choose CellBase (behavior data stored in CellBase structure)
choosecb(cbname)
allsessions = listtag('allsessions'); % list all sessions
NumSessions = size(allsessions,1); % number of all sessions of all animals

if ~isfolder(resdir) % if results directory does not exist, make it
    mkdir(resdir)
end

if isbehav % Need for behavior data preprocessing
for n = 1:NumSessions
    ifsave=1;
    filepath = fullfile(getpref('cellbase', 'datapath'), allsessions{n,1}, allsessions{n,2}, [allsessions{n,1} allsessions{n,2} '.mat']);
    TE = solo2trialevents_pavlovian_FF(filepath,ifsave)
end
end

% Time window
wn = [-3 3];   % in seconds
twin = [0 1]; % test window for statistics - 1s after cue onset
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector
tinx = [find(time==twin(1)) find(time==twin(2))]; % find test window indices in time vector

% PETH for likely reward and punishment
[Hit_allpsth, FA_allpsth, Hit_spt, FA_spt] = deal([]);
for iS = 1:NumSessions
    sessionid = allsessions(iS,:);
    [spsth_hit, spt_hit, spsth_fa, spt_fa] = main(sessionid,wn,dt);
    Hit_allpsth = [Hit_allpsth; spsth_hit]; % Likely reward trials PETH
    FA_allpsth = [FA_allpsth; spsth_fa]; % Likely punishment trials PETH
    Hit_spt = [Hit_spt; mean(sum(spt_hit(:,tinx(1):tinx(2)),2))]; % Likely reward trials lick raster
    FA_spt = [FA_spt; mean(sum(spt_fa(:,tinx(1):tinx(2)),2))]; % Likely punishment trials lick rasters
end

% Boxstat to compare lick rate accross likely reward and likely punishment trials
boxstat2(Hit_spt, FA_spt, 'T1', 'T2', 0.05, 'paired')
fnm = fullfile(resdir,'Boxstat_lickPSTH.fig');
saveas(gcf,fnm)
fnm = fullfile(resdir,'Boxstat_lickPSTH.jpeg');
saveas(gcf,fnm)

% Bar plot
figure;
bar([1 2], [median(Hit_spt) median(FA_spt)])
sem1 = se_of_median(Hit_spt);
sem2=se_of_median(FA_spt);
hold on
errorbar([1 2], [median(Hit_spt) median(FA_spt)], [sem1 sem2], [sem1 sem2])
fnm = fullfile(resdir,'Barplot_lickPSTH.fig');
saveas(gcf,fnm)
set(gcf,'renderer','painters')
fnm = fullfile(resdir,'Barplot_lickPSTH.eps');
saveas(gcf,fnm)

% PETH Plot & save
H = figure;
green = [51 204 51] / 255;   % colors for plotting
red = [216 41 0] / 255;
errorshade(time,nanmean(Hit_allpsth),nanstd(Hit_allpsth)/sqrt(size(Hit_allpsth,1)),...
    'LineColor',green,'ShadeColor',green)
hold on
errorshade(time,nanmean(FA_allpsth),nanstd(FA_allpsth)/sqrt(size(FA_allpsth,1)),...
    'LineColor',red,'ShadeColor',red)
set(H,'renderer', 'painters')
fnm = fullfile(resdir,'average_lickPSTH.fig');
saveas(H,fnm)
fnm = fullfile(resdir,'average_lickPSTH.eps');
saveas(H,fnm)
keyboard

% -------------------------------------------------------------------------
function [spsth_hit, spt_hit, spsth_fa, spt_fa] = main(cellid,wn,dt)

% Filter input
filterinput_hit = 'TrialType==1';
filterinput_fa = 'TrialType==2';

% Calcualte lick PSTH
[~, spsth_hit, ~, ~, spt_hit] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.05,'event_filter','custom','filterinput',filterinput_hit,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_fa, ~, ~, spt_fa] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.05,'event_filter','custom','filterinput',filterinput_fa,...
    'isadaptive',0,'maxtrialno',Inf);

% Lick raster
H = figure;
viewlick(cellid,'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','sigma', 0.05,...
    'ShowEvents',{{'StimulusOn' 'StimulusOff'}},...
    'Partitions','#TrialType','window',wn)
maximize_figure(H)
close(H)