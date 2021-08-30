function Ach_lick_poppsth(achcells,resdir,issave)
%ACH_LICK_POPPSTH(ACHCELLS,RESDIR,ISSAVE) color coded lick PETH.
%   ACH_LICK_POPPSTH plots color coded lick PETH for each recording
%   session where cells defined in ACHCELLS were recorded. PETHs are aligned
%   to cue onset. Plots are saved to RESDIR if ISSAVE is true.

%   See also ULTIMATE_PSTH

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   29-Apr-2020

%   Code review: BH 5/11/20

% Input arguments
if nargin < 1 || isempty(achcells) % if no cellids defined in the input
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

% Sessions
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
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector

[Hit_allpsth, FA_allpsth] = deal([]); % preallocate space for lick psths respectively for reward and punishment cue
for iS = 1:NumSessions
    cSession = Sessions{iS}; % extract sessionIDs in the form of 1x2 cell
    [rat,remain] = strtok(cSession,'_');
    [session, ~] = strtok(remain, '_');
    sessionid = [{rat}, {session}];
    [spsth_hit, spsth_fa] = main(sessionid,wn,dt); % calculate PSTHs
    Hit_allpsth = [Hit_allpsth; spsth_hit];
    FA_allpsth = [FA_allpsth; spsth_fa];
end

% Color coded lick raster
PSTH_all = [{Hit_allpsth} {FA_allpsth}];
cax = [0 0]; % inital values for colorbar
for j = 1:length(PSTH_all)
    if j==1  % sort based on T1 response latency
    [m1, m2] = max(PSTH_all{j},[],2);
    [srt, Ia] = sort(m1,'descend');  
    end
    subplot(1,2,j);
    imagesc(time,1:size(PSTH_all{j},1),PSTH_all{j}(Ia,:)) % plot
    colormap(jet)
    colorbar
    c = caxis;  % use same colorbar scale for subplots
    c1 = [min([c, cax]), max([c, cax])];
    caxis(c1);
    cax = c1;
end

% Save figure (.fig, .jpeg and .eps format)
if issave 
    set(gcf, 'Renderer', 'painters');
    fnm = fullfile(resdir, 'Lick_poppsth.fig');
    saveas(gcf,fnm)
    fnm2 = fullfile(resdir, 'Lick_poppsth.jpg');
    saveas(gcf,fnm2)
    fnm3 = fullfile(resdir, 'Lick_poppsth.eps');
    saveas(gcf,fnm3)   
end

% -------------------------------------------------------------------------
function [spsth_hit, spsth_fa] = main(cellid,wn,dt)

% Filter input
filterinput_hit = 'TrialType==1';
filterinput_fa = 'TrialType==2';

% Calcualte lick PSTH
[~, spsth_hit] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.02,'event_filter','custom','filterinput',filterinput_hit,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_fa] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.02,'event_filter','custom','filterinput',filterinput_fa,...
    'isadaptive',0,'maxtrialno',Inf);