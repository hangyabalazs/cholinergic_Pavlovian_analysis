function Ach_lick_psth_summary(achcells,resdir,issave)
%ACH_LICK_PSTH_SUMMARY(ACHCELLS, RESDIR, ISSAVE)   Average lick PETH.
%   ACH_LICK_PSTH_SUMMARY plots average lick PETH (beam break time stamps
%   aligned to an event). All sessions from mice are used.
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   29-Apr-2020

%   Code review: BH 5/11/20

% Input arguments
if nargin < 1 || isempty(achcells)
    achcells = select_ach_cells(achcells);
end

if nargin < 2
    resdir = ('D:\pavlovian_cholinergic_cellbase\_paper_figs\Fig1_behavior\');   % result directory
end
if nargin < 3
    issave = true;
end

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Select sessions: mice w no feedback delay that had enough training - last 10
% sessions (only first sessions if there were more sessions per day)
cellids = achcells'; % BF cholinergic cell IDs
NumCells = length(cellids); % number of cells
Sessions =cell(1,NumCells); % sessions of recording
for i = 1:NumCells % extract session names from cellids
    [ratname,session,~,~] = cellid2tags(cellids{i});
    Sessions{i} = [ratname '_' session];
end
Sessions = unique(Sessions); % eliminate duplications
NumSessions = length(Sessions);

% Time window
wn = [-3 3];   % in seconds
twin = [0 1];
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector
tinx = [find(time==twin(1)) find(time==twin(2))];

% PETH
[Hit_allpsth, FA_allpsth, Hit_allspt, FA_allspt] = deal([]);
for iS = 1:NumSessions
    cSession = Sessions{iS}; % extract sessionIDs in the form of 1x2 cell
    [rat,remain] = strtok(cSession,'_');
    [session, ~] = strtok(remain, '_');
    sessionid = [{rat}, {session}];
    [spsth_hit, spt_hit, spsth_fa spt_fa] = main(sessionid,wn,dt);
    Hit_allpsth = [Hit_allpsth; spsth_hit];
    FA_allpsth = [FA_allpsth; spsth_fa];
    Hit_allspt = [Hit_allspt; mean(sum(spt_hit(:,tinx(1):tinx(2)),2))];
    FA_allspt = [FA_allspt; mean(sum(spt_fa(:,tinx(1):tinx(2)),2))];
    
    H = gcf;
    if issave
        cellidt = [sessionid{1} '_' sessionid{2}];
        fnm = fullfile(resdir, [cellidt '_LICK.fig']);   % save
        saveas(H,fnm)
        fnm = fullfile(resdir, [cellidt '_LICK.jpg']);
        saveas(H,fnm)
    end
    close(H)
end

boxstat2(Hit_allspt, FA_allspt, 'T1', 'T2', 0.05, 'paired')
bar([1 2], [median(Hit_allspt) median(FA_allspt)])
sem1 = se_of_median(Hit_allspt);
sem2=se_of_median(FA_allspt);
hold on
errorbar([1 2], [median(Hit_allspt) median(FA_allspt)], [sem1 sem2], [sem1 sem2])

% Plot & save
H = figure;
green = [51 204 51] / 255;   % colors for plotting
red = [216 41 0] / 255;
errorshade(time,nanmean(Hit_allpsth),nanstd(Hit_allpsth)/sqrt(size(Hit_allpsth,1)),...
    'LineColor',green,'ShadeColor',green)
hold on
errorshade(time,nanmean(FA_allpsth),nanstd(FA_allpsth)/sqrt(size(FA_allpsth,1)),...
    'LineColor',red,'ShadeColor',red)
if issave
    set(H,'renderer', 'painters')
    fnm = fullfile(resdir,'average_lickPSTH.fig');
    saveas(H,fnm)
    fnm = fullfile(resdir,'average_lickPSTH.eps');
    saveas(H,fnm)
end

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