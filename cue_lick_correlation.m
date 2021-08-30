function cue_lick_correlation(cellids, resdir)
%CUE_LICK_CORRELATION   Lick PSTH partitioned by cue response.
%   CUE_LICK_CORRELATION(CELLIDS, RESDIR) compares anticipatory licking
%   when the cue response was high or low. Behvaioral sessions of cells
%   (CELLIDS) are used. Results are saved in RESDIR. Lick PSTH partitioned
%   by cue response (separately for the two cues) and reaction time
%   comparisons for cue response median split are generated.
%
%   See also AVG_PSTH_CHOLINERGIC and ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine
%   hegedus.panna@koki.mta.hu

%   Code review: BH 11/26/20

% Directories
if ~isfolder(resdir) % if results directory does not exist
    mkdir(resdir)
end

% Time variables
wn = [-4 4]; % time window
dt = 0.001; % time resolution
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector

twin = [0 0.5];  % test window for statistics
bwin = [-1 0];   % baseline window

% Index for time 0
st = abs(wn(1)) / dt;   % in ms
nullindex = st + 1;

% PSTHs for cue response (to generate bin raster)
R1 = runanalysis(@ultimate_psth,...
    'trial', 'StimulusOn', wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',0,'forcesmoothedstat',true,...
    'maxtrialno',Inf,'baselinewin',bwin,'testwin',twin,'relative_threshold',0.1,'cellids',cellids);

% Window for testing the potential effect
WNt = [twin(1)/dt+nullindex twin(2)/dt+nullindex-1];   % test window; convert to indices
WNt = round(WNt);
lWNt = WNt(2) - WNt(1) + 1;   % length of test window
BR1 = R1(:,5); % binraster for cue1
NumCell = length(cellids);
[RT1_T1, RT2_T1, RT1_T2, RT2_T2] = deal(cell(1,NumCell));

for iC = 1:NumCell % count spikes for every session for every cell, do median split and modify TrialEvents struct
    
    binrast = BR1{iC};  % bin raster
    
    avg_spikecount_cell = nan(1,size(binrast ,1)); % calculate average spike counts from binraster
    for j = 1:size(binrast ,1)   % loop through trials
        avg_spikecount_cell(j) = sum(binrast(j, WNt(1):WNt(2)), 2) / lWNt / dt;
    end
    
    TE = loadcb(cellids{iC}, 'TrialEvents'); % load TrialEvents struct
    T1_inx = (TE.TrialType==1);
    T2_inx = (TE.TrialType==2);
    
    m1 = median(avg_spikecount_cell(T1_inx)); % calculate median for both Trialtype 1 and 2
    m2 = median(avg_spikecount_cell(T2_inx));
    
    TE.Cueresp = nan(1,length(TE.NTrials)); % make new variable in TE struct for cue response partition
    TE.Cueresp(T1_inx & (avg_spikecount_cell>m1)) = 1; % if bigger than median
    TE.Cueresp(T1_inx & (avg_spikecount_cell<=m1)) = 2; % if smaller than median
    TE.Cueresp(T2_inx & (avg_spikecount_cell>m2)) = 1;
    TE.Cueresp(T2_inx & (avg_spikecount_cell<=m2)) = 2;
    RT1_T1{iC} = TE.RT(T1_inx & avg_spikecount_cell>m1); % reaction times grouped by median split of cue response (separately for T1 and T2)
    RT2_T1{iC} = TE.RT(T1_inx & avg_spikecount_cell<=m1);
    RT1_T2{iC} = TE.RT(T2_inx & avg_spikecount_cell>m2);
    RT2_T2{iC} = TE.RT(T2_inx & avg_spikecount_cell<=m2);
    
    aID = getvalue('RatID_tag', cellids{iC});
    sID = getvalue('SessionID_tag', cellids{iC});
    savedir = fullfile(getpref('cellbase', 'datapath'), aID{1}, sID{1});
    savename = 'TrialEvents.mat'; % save TE struct
    savename2 = 'TE.mat'; % save TE struct
    save(fullfile(savedir, savename),'-struct','TE');
    save(fullfile(savedir, savename2),'-struct','TE');
end

% Bar graph showing RT partitioned to big and small cue response
RT1_T1_m = cellfun(@nanmean,RT1_T1);  % average within cells
RT2_T1_m = cellfun(@nanmean,RT2_T1);
RT1_T2_m = cellfun(@nanmean,RT1_T2);
RT2_T2_m = cellfun(@nanmean,RT2_T2);

% Statistics
[~, p1] = boxstat(RT1_T1_m,RT2_T1_m,'High cue resp','Low cue resp',0.05,'paired');
filename = fullfile(resdir,'box_RT_T1_cueresp.jpg');
saveas(gcf,filename);
set(gcf, 'renderer', 'painters');
filename = fullfile(resdir,'box_RT_T1_cueresp.eps');
saveas(gcf,filename);
close(gcf)
[~, p2] = boxstat(RT1_T2_m,RT2_T2_m,'High cue resp','Low cue resp',0.05,'paired');
filename = fullfile(resdir,'box_RT_T2_cueresp.fig');
saveas(gcf,filename);
set(gcf, 'renderer', 'painters');
filename = fullfile(resdir,'box_RT_T2_cueresp.eps');
saveas(gcf,filename);
close(gcf)

% Bar graph and SE for RT after reward predicting cue
B1 = figure;
data1 = [nanmedian(RT1_T1_m) nanmedian(RT2_T1_m)];
bar([1,2], data1, 'FaceColor', [0.9 0.9 0.9]);
hold on
bootstrap_se1 = se_of_median(RT1_T1_m);
bootstrap_se2 = se_of_median(RT2_T1_m);
er = errorbar([1,2],data1,[bootstrap_se1 bootstrap_se2],'k.');
line([1 2], [RT1_T1_m', RT2_T1_m'],'Color', [0.7 0.7 0.7])

% Save figure
set(B1, 'Renderer', 'painters');
filename = fullfile(resdir,'compare_RT_T1_cueresp.jpg');
saveas(B1,filename);
filename = fullfile(resdir,'compare_RT_T1_cueresp.eps');
saveas(B1,filename);
close(B1)

% Bar graph and SE for RT after punishment predicting cue
B2 = figure;
data2 = [nanmedian(RT1_T2_m) nanmedian(RT2_T2_m)];
bar([1,2], data2, 'FaceColor', [0.9 0.9 0.9]);
hold on
bootstrap2_se1 = se_of_median(RT1_T2_m);
bootstrap2_se2 = se_of_median(RT2_T2_m);
er = errorbar([1,2],data2,[bootstrap2_se1 bootstrap2_se2],'k.');
line([1 2], [RT1_T2_m', RT2_T2_m'],'Color', [0.7 0.7 0.7])

% Save figure
set(B2, 'Renderer', 'painters');
filename = fullfile(resdir,'compare_RT_T2_cueresp.fig');
saveas(B2,filename);
filename = fullfile(resdir,'compare_RT_T2_cueresp.eps');
saveas(B2,filename);
close(B2)

% Convert cellids to sessionids - for lick psth
animalIDs = getvalue('RatID_tag', cellids);
sessionID = getvalue('SessionID_tag', cellids);
sessionids = cell(1,NumCell);
for k = 1:NumCell
    sessionids{k} = {animalIDs{k} sessionID{k}};
end

% Average PSTHs partitioned to big and small cue response
avg_psth_cholinergic(sessionids, 'cue', 'cue_corr', resdir) % make average PSTH