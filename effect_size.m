function effect_size(achcells, responsetype, ispart, resdir, data_type)
%EFFECT_SIZE quantify effect size.
%   
%   EFFECT_SIZE(ACHCELLS, RESPONSETYPE, ISPART, RESDIR, DATA_TYPE)compare
%   cue, reward and punishment (RESPONSETYPE) response size of ACHCELLS 
%   under expected and unexpected conditions. If ISPART is true, trials are 
%   partitioned by the defined variable and plotted separately. Results are
%   saved to RESDIR.

%   EFFECT_SIZE(ACHCELLS, RESPONSETYPE, ISPART, RESDIR, DATA_TYPE)
%   calls ULTIMATE_PSTH for virtual (simulated) spikes if DATA_TYPE input
%   argument is 'virtual'. This can be used to analyze simulated data the
%   same way as real spiking.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   01-May-2020

% Input arguments
if nargin < 1 || isempty(achcells)
    achcells = select_ach_cells(achcells);
end
NumCell = length(achcells);
if nargin < 5
    data_type = 'real';  % control simulated data
end

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Define filterinput
switch responsetype
    case 'cue' % cure response
        filters = {'TrialType==1' 'TrialType==2'};
        align = 'StimulusOn';
        twin = [0 0.5];
        bwin = [-1 0];   % baseline window
    case 'lickcue' % cue reaponse with anticipatory licking
        filters = {'Lick==1&TrialType==1' 'Lick==1&TrialType==2'};
        align = 'StimulusOn';
        twin = [0 0.5];
        bwin = [-1 0];   % baseline window
    case 'reward' % reward response
        filters = {'AllReward==1' 'AllReward==2'};
        align = 'DeliverAllFeedback';
        twin = [0 0.2];
        bwin = [-3 -2];   % baseline window
    case 'lickreward' % reward response after anticipatory licking
        filters = {'Reward==1' 'Reward==2'};
        align = 'DeliverAllFeedback';
        twin = [0 0.2];
        bwin = [-3 -2];   % baseline window
    case 'lickpunish' % lick & punish response
        filters = {'FalseAlarm==1&TrialType==1' 'FalseAlarm==1&TrialType==2'};
        align = 'DeliverAllFeedback';
        twin = [0 0.2];
        bwin = [-3 -2];   % baseline window
    case 'punish' % punishment response
        filters = {'Punishment==1' 'Punishment==2'};
        align = 'DeliverAllFeedback';
        twin = [0 0.2];
        bwin = [-3 -2];   % baseline window
    case 'omission'  % omission responses
        filters = {'Omission==1' 'Omission==2'};  % omission after likely reward (1) or likely punishment cue (2)
        align = 'DeliverAllFeedback';
        twin = [0 0.5];
        bwin = [-3 -2];   % baseline window
end

% Filter trials for anticipatory lick response
switch ispart
    case 'lick'
        part = '#Lick';
        evtype = 'trial';
    case 'cue_corr'
        part = '#Cueresp';
        evtype = 'lick';
        twin = [0 1];  % licks can be counted in the entire period before reinforcement
        bwin = [-2 0];   % baseline window
    case 'cue'
        part = '#Cueresp';
        evtype = 'trial';
    case 'previous_outcome'
        part = '#PreviousOutcome';
        evtype = 'trial';
    case 'none'
        part = 'all';
        evtype = 'trial';
end

% Time window
wn = [-4 4];
dt = 0.001;
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector

% Index for time 0
st = abs(wn(1)) / dt;   % in ms
nullindex = st + 1;

% PSTH
R1 = runanalysis(@ultimate_psth,...
    evtype, align, wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts',part,'isadaptive',0,...
    'event_filter','custom', 'filterinput',filters{1},'maxtrialno',Inf,'data_type',data_type,...
    'baselinewin',bwin,'testwin',twin, 'forcesmoothedstat', true, 'relative_threshold',0.1,'cellids',achcells);
R2 = runanalysis(@ultimate_psth,...
    evtype, align, wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts',part,'isadaptive',0,...
    'event_filter','custom', 'filterinput',filters{2},'maxtrialno',Inf,'data_type',data_type,...
    'baselinewin',bwin,'testwin',twin, 'forcesmoothedstat', true, 'relative_threshold',0.1,'cellids',achcells);

% Window for testing the potential effect
WNt = [twin(1)/dt+nullindex twin(2)/dt+nullindex-1];   % test window; convert to indices
WNt = round(WNt);
lWNt = WNt(2) - WNt(1) + 1;   % length of test window

BR1 = R1(:,5); % binraster for cue1
BR2 = R2(:,5); % binraster for cue2

    [psth_R1, psth_R2] = deal(nan(NumCell,size(R1{1,1},2)));
    baseline_inx = find(time==bwin(1)*1000):find(time==bwin(2)*1000); % baseline indices
    Rinx = 2;
    for iC = 1:NumCell % Z-score by baseline
    
            T1size = size(R1{iC,5},1);
            T2size = size(R2{iC,5},1);
     
        if (T1size>5) && (T2size>5) %  only use cells with >5 trials in each partition
            
            mn = mean(R1{iC,Rinx}(1,baseline_inx)); % means and sds calculated for the average PSTH
            sd = std(R1{iC,Rinx}(1,baseline_inx));
            
                mn2 = mn;
                sd2 = sd;

           
                psth_R1(iC,:) = (R1{iC,Rinx} - mn) / sd;   % the same mean and SD used within a cell,
                psth_R2(iC,:) = (R2{iC,Rinx} - mn2) / sd2;
            
        end        
    end
    psth_R1(any(isnan(psth_R1), 2), :) = []; % remove NaNs
    psth_R2(any(isnan(psth_R2), 2), :) = [];
    
    psth_R1_t=max(mean(psth_R1(:,[WNt(1):WNt(2)])));
    psth_R2_t=max(mean(psth_R2(:,[WNt(1):WNt(2)])));

diff_all=(psth_R1_t./psth_R2_t)*100-100;
diff_all2=(psth_R2_t./psth_R1_t)*100-100;

% Save results
save(fullfile(resdir, [responsetype '_effect_size.mat']), 'diff_all', 'diff_all2');