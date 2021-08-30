function compare_cue_RT(cellids, resdir)
%COMPARE_CUE_RT   Cue response latency and reacion time comparison.
%   COMPARE_CUE_RT(CELLIDS, RESDIR) compares average cue response of cells
%   listed in CELLIDS and corresponding average reaction time. PSTHs are
%   saved to RESDIR.
%
%   See also VIEWLICK, VIEWCELL2B and ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine
%   hegedus.panna@koki.mta.hu

%   Code review: BH 7/30/20, 12/1/20

% Results directory
if ~isfolder(resdir) % make results directory
    mkdir(resdir)
end

% Cell IDs
if isempty(cellids) % if no input cellids are given, select all cholinergic neurons
    cellids = select_ach_cells(cellids);
end

% Parameters
wn = [-1 2.5]; % Time window for PSTH
bwin = [-1 0]; % baseline window
dt = 0.001; % time resolution
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector
sigma = 0.08;

% Index for time 0
st = abs(wn(1)) / dt;   % in ms
nullindex = st + 1;
NumCell = length(cellids);

% Call for single example or average
if length(cellids) == 1 % if examine only one cell
    
    % Lick times
    TE = loadcb(cellids,'TrialEvents');
    margin = sigma * 3;     % add an extra margin to the windows
    time = wn(1)-margin:dt:wn(2)+margin;   % time base array
    NumTrials = length(TE.LickIn); % number of trials
    stimes = arrayfun(@(s)TE.LickIn{s}-TE.StimulusOn(s),1:NumTrials,'UniformOutput',false); % stimulus times
%     if ~isfield(TE,'FirstLickAfterStimulus')  % add lick variables to TrialEvents
        lickvariables(cellids,TE,NumTrials)
%     end
    
    % Spike raster
    pause(0.01)
    viewcell2b(cellids,'TriggerName','StimulusOn','SortEvent','FirstLickAfterStimulus','sigma', sigma,...
        'eventtype','behav','ShowEvents',{{'FirstLickAfterStimulus'}},'Partitions','#GoLick','window',wn)
    
    cellidt = cellids{1}; % save psth
    cellidt(cellidt=='.') = '_';
    set(gcf, 'renderer', 'painters');
    fnm = fullfile(resdir,[cellidt '_StimulusOn_example.eps']);   % save
    fnm2 = fullfile(resdir,[cellidt '_StimulusOn_example.jpg']);   % save
    fnm3 = fullfile(resdir,[cellidt '_StimulusOn_example.fig']);   % save
    saveas(gcf,fnm)
    saveas(gcf,fnm2)
    saveas(gcf,fnm3)
    close(gcf)
    
    % Lick raster
    H = figure;
    aID = getvalue('RatID_tag',cellids{1});
    sID = getvalue('SessionID_tag',cellids{1});
    viewlick({aID{1}, sID{1}},'TriggerName','StimulusOn','SortEvent','FirstLickAfterStimulus','sigma', sigma,...
        'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#GoLick','window',wn);
    set(H, 'renderer', 'painters');
    
    fnm = fullfile(resdir,[cellidt '_LickRaster_example.eps']);   % save
    fnm2 = fullfile(resdir,[cellidt '_LickRaster_example.jpg']);   % save
    fnm3 = fullfile(resdir,[cellidt '_LickRaster_example.fig']);   % save
    saveas(H,fnm)
    saveas(H,fnm2)
    saveas(H,fnm3)
    close(H)
    
else % if there is a cell array of cellids
    
    % PSTH for cue response
    C1 = runanalysis(@ultimate_psth,...
        'trial', 'StimulusOn', wn,...
        'dt',0.001,'display',false,'sigma',sigma,'isadaptive',0,...
        'event_filter','custom', 'filterinput','TrialType==1','maxtrialno',Inf,...
        'baselinewin',bwin,'testwin',[0 0.5],'relative_threshold',0.1,'cellids',cellids);
    
    C2 = runanalysis(@ultimate_psth,...
        'trial', 'StimulusOn', wn,...
        'dt',0.001,'display',false,'sigma',sigma,'isadaptive',0,...
        'event_filter','custom', 'filterinput','TrialType==2','maxtrialno',Inf,...
        'baselinewin',[-1 0],'testwin',[0 0.5],'relative_threshold',0.1,'cellids',cellids);
    
    [psth_C1, psth_C2] = deal(nan(NumCell,size(C1{1,1},2)));
    baseline_inx = find(time==bwin(1)*1000):find(time==bwin(2)*1000); % baseline indices
    Rinx = 2;
    
    for iC = 1:NumCell % Z-score by baseline
        if size(C1{iC,5},1) > 5 && size(C2{iC,5},1) > 5 %  only use cells with >5 trials in each partition
            mn = mean(C1{iC,Rinx}(1,baseline_inx)); % means and sds calculated for the average PSTH
            sd = std(C1{iC,Rinx}(1,baseline_inx));
            psth_C1(iC,:) = (C1{iC,Rinx} - mn) / sd;   % the same mean and SD used within a cell,
            psth_C2(iC,:) = (C2{iC,Rinx} - mn) / sd;
        end
    end
    psth_C1(any(isnan(psth_C1), 2), :) = []; % remove NaNs
    psth_C2(any(isnan(psth_C2), 2), :) = [];
    
    % Convert cellids to sessionids - needed for lick psth
    animalIDs = getvalue('RatID_tag', cellids);
    sessionID = getvalue('SessionID_tag', cellids);
    
    sessionids = cell(1,NumCell);
    for k = 1:NumCell
        sessionids{k} = {animalIDs{k} sessionID{k}};
    end
    
    % Lick PSTH
    L1 = runanalysis(@ultimate_psth,...
        'lick','StimulusOn',wn,...
        'dt',dt,'sigma',sigma,'event_filter','custom','filterinput','TrialType==1',...
        'isadaptive',0,'maxtrialno',Inf, 'cellids', sessionids);
    
    L2 = runanalysis(@ultimate_psth,...
        'lick','StimulusOn',wn,...
        'dt',dt,'sigma',sigma,'event_filter','custom','filterinput','TrialType==2',...
        'isadaptive',0,'maxtrialno',Inf, 'cellids',sessionids);
    
    [psth_L1, psth_L2] = deal(nan(NumCell,size(L1{1,1},2)));
    
    for iC = 1:NumCell
        
        % No normalization for lick PETH for consistency
        if size(L1{iC,5},1) > 5 && size(L2{iC,5},1) > 5 %  only use cells with >5 trials in each partition
            psth_L1(iC,:) = L1{iC,Rinx};
            psth_L2(iC,:) = L2{iC,Rinx};
        end
    end
    psth_L1(any(isnan(psth_L1), 2), :) = []; % remove NaNs
    psth_L2(any(isnan(psth_L2), 2), :) = [];
    
    % Plot & save
    H = figure;
    subplot(2,1,1);
    green = [51 204 51] / 255;   % color for neuronal response plotting
    yellow = [229 245 16] / 255; % color for lick response plotting
    errorshade(time,nanmean(psth_C1),nanstd(psth_C1)/sqrt(size(psth_C1,1)),...
        'LineColor',green,'ShadeColor',green)
    subplot(2,1,2);
    errorshade(time,nanmean(psth_L1),nanstd(psth_L1)/sqrt(size(psth_L1,1)),...
        'LineColor',yellow,'ShadeColor',yellow)
    set(H, 'renderer', 'painters');
    fnm = fullfile(resdir,'TrialType1_LickRasterand_PSTH_all.eps');   % save
    fnm2 = fullfile(resdir,'TrialType1_LickRasterand_PSTH_all.jpg');   % save
    fnm3 = fullfile(resdir,'TrialType1_LickRasterand_PSTH_all.fig');   % save
    
    saveas(H,fnm)
    saveas(H,fnm2)
    saveas(H,fnm3)
    close(H)
    
    % Plot & save
    G = figure;
    subplot(2,1,1);
    red = [216 41 0] / 255;
    errorshade(time,nanmean(psth_C2),nanstd(psth_C2)/sqrt(size(psth_C2,1)),...
        'LineColor',red,'ShadeColor',red)
    subplot(2,1,2);
    errorshade(time,nanmean(psth_L2),nanstd(psth_L2)/sqrt(size(psth_L2,1)),...
        'LineColor',yellow,'ShadeColor',yellow)
    set(G, 'renderer', 'painters');
    fnm = fullfile(resdir,'TrialType2_LickRasterand_PSTH_all.eps');   % save
    fnm2 = fullfile(resdir,'TrialType2_LickRasterand_PSTH_all.jpg');   % save
    fnm3 = fullfile(resdir,'TrialType2_LickRasterand_PSTH_all.fig');   % save
    
    saveas(G,fnm)
    saveas(G,fnm2)
    saveas(G,fnm3)
    close(G)
end

% -------------------------------------------------------------------------
function lickvariables(cellid,TE,NumTrials)

TE.FirstLickAfterStimulus = nan(1,NumTrials);
for iT = 1:NumTrials
    if ~isempty(TE.LickIn{iT})   % only if there were licks after stimulus
        TE.FirstLickAfterStimulus(1,iT) = TE.LickIn{iT}(...
            find(TE.LickIn{iT}>TE.StimulusOn(iT),1,'first'));
    end
end
TE.FirstLick = TE.FirstLickAfterStimulus;
[TE.GoLick, TE.NoGoLick, TE.Lick] = deal(nan(1,NumTrials));
TE.GoLick(~isnan(TE.GoRT)) = 1;  % lick response in TrialType1 trials
TE.NoGoLick(~isnan(TE.NoGoRT)) = 1;  % lick response in TrialType2 trials

TE.Lick(~isnan(TE.GoRT))=1;
TE.Lick(~isnan(TE.NoGoRT))=1;
TE.Lick(isnan(TE.GoRT)&isnan(TE.NoGoRT))=2;

aID = getvalue('RatID_tag', cellid);
sID = getvalue('SessionID_tag', cellid);
savedir = fullfile(getpref('cellbase', 'datapath'), aID{1}, sID{1});
savename = 'TrialEvents.mat'; % save TE struct
savename2 = 'TE.mat'; % save TE struct
save(fullfile(savedir, savename),'-struct','TE');
save(fullfile(savedir, savename2),'-struct','TE');