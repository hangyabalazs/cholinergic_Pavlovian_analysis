function cholinergic_LS(cellids, resdir)
%CHOLINERGIS_LS Light-alinged PSTH.
%   CHOLINERGIC_LS(CELLIDS, RESDIR) Calculates light-aligned PSTH with
%   small, [-0.005 0.005] s window for neurons listed in CELLIDS. The PSTHs
%   are saved to RESDIR.

%   See also MAKESTIMEVENTS2_P, PREALIGNSPIKES, VIEWCELL2b and
%   QUICKANALYSIS_PAVLOVIAN2_P

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   16-Apr-2020

%   Code review: BH 5/11/20

% Input check
narginchk(0,3)
if nargin < 1
    cellids = [];
end
if nargin < 2
    resdir = fullfile(getpref('cellbase', 'datapath'), 'Cholinergic_LS\'); %Results directory
end

% Results directory
if ~isfolder(resdir)
    mkdir(resdir)
end

% Create stimulus events
for i = 1:length(cellids)
    currentCell = cellids{i};
    animalID = currentCell(1:5);
    sessionID = currentCell(7:13);
    fullpth = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
    
    MakeStimEvents2_p(fullpth,'BurstStartNttl',4)
end

% Prealign spikes to stimulus events
problem_stim_cellid = [];
for iC = 1:length(cellids)
    cellid = cellids(iC);
    try
        prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim,'filetype','stim','ifsave',1,'ifappend',0)
    catch
        disp('Error in prealignSpikes.');
        problem_stim_cellid = [problem_stim_cellid cellid];
    end
end

% View light-triggered raster and PSTH
TrigEvent = 'BurstOn'; % Inputs for viewcell2b (for detailed information, see VIEWCELL2B.M)
SEvent = 'BurstOff';
win = [-0.005 0.005];
parts = '#BurstNPulse';
dt = 0.001;
sigma = 0.001;
PSTHstd = 'on';
ShEvent = {{'PulseOn','PulseOff','BurstOff'}};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
for iCell = 1:length(cellids)
    cellid = cellids(iCell);
    figure;
    viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','off')
    
    cellidt = cellid{1};
    cellidt(cellidt=='.') = '_';
    set(gcf, 'Renderer', 'painters')
    fnm = fullfile(resdir,[cellidt '_LS.fig']);   % save
    saveas(gcf,fnm)
    fnm = fullfile(resdir,[cellidt '_LS.jpg']);   % save
    saveas(gcf,fnm)
    close(gcf)
end