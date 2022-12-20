function cueresp_RT_regression(cellids, resdir)
% CUERESP_RT_REGRESSION Linear regression of cue responses, reaction times
% and anticipatory lick rates
%
% CUERESP_RT_REGRESSION(CELLIDS, RESDIR) Performs linear regression of cue
% responses, anticipatory lick rates and reaction times.

% See also POLYPREDCICALL

% Panna Hegedus
% Institute of Experimental Medicine
% hegedus.panna@koki.mta.hu

if ~isfolder(resdir)
    mkdir(resdir);
end

NumCells = length(cellids); % Number of cells
[RT1_avg, RT2_avg, LR1_avg, LR2_avg, FR1_avg, FR2_avg] = deal(nan(1,NumCells)); % Preallocate space for avg lick rate and reaction time

align = 'StimulusOn';
twin = [0 0.5];
bwin = [-1 0];   % baseline window

% Time window
wn = [-4 4];
dt = 0.001;
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector

% Index for time 0
st = abs(wn(1)) / dt;   % in ms
nullindex = st + 1;

% Window for testing the potential effect
WNt = [twin(1)/dt+nullindex twin(2)/dt+nullindex-1];   % test window; convert to indices
WNt = round(WNt);
lWNt = WNt(2) - WNt(1) + 1;   % length of test window


R1 = runanalysis(@ultimate_psth,...
    'trial', 'StimulusOn', wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',0,'forcesmoothedstat',true,...
    'maxtrialno',Inf,'baselinewin',bwin,'testwin',twin,'relative_threshold',0.1,'cellids',cellids);
BR = R1(:,5); % binraster for cue

for iC = 1:NumCells
    % Calculate reaction time and lick rate
    TE = loadcb(cellids{iC}, 'TrialEvents'); % load TrialEvents struct
    T1_inx = (TE.TrialType==1); % indices of T1 trials
    T2_inx = (TE.TrialType==2);  % indices of T2 trials
    RT1_avg(iC) = nanmean(TE.RT(T1_inx));
    RT2_avg(iC) = nanmean(TE.RT(T2_inx));
    LR1_avg(iC) = calc_avg_LR(TE.LickIn(T1_inx), TE.StimulusOn(T1_inx));
    LR2_avg(iC) = calc_avg_LR(TE.LickIn(T2_inx), TE.StimulusOn(T2_inx));
    
    % Calculate cue response
    binrast = BR{iC};  % bin raster
    
    avg_spikecount_cell = nan(1,size(binrast,1)); % calculate average spike counts from binraster
    for j = 1:size(binrast ,1)   % loop through trials
        avg_spikecount_cell(j) = sum(binrast(j, WNt(1):WNt(2)), 2) / lWNt / dt;
    end
    
    FR1_avg(iC) = nanmean(avg_spikecount_cell(T1_inx)); % calculate mean FR for both Trialtype 1 and 2
    FR2_avg(iC) = nanmean(avg_spikecount_cell(T2_inx));
end

% Calculate regression
% Regression of lick rate and reaction time
polypredcicall(FR1_avg,RT1_avg,0.95,'robust')
xlabel('FR to reward cue');
ylabel('RT to reward cue');
fnm = fullfile(resdir,'FR_vs_RT_1.fig'); % save figures
saveas(gcf,fnm);
fnm = fullfile(resdir,'FR_vs_RT_1.jpg');
saveas(gcf,fnm);
set(gcf, 'renderer', 'painters')
fnm = fullfile(resdir,'FR_vs_RT_1.eps');
saveas(gcf,fnm);

polypredcicall(FR2_avg,RT2_avg,0.95,'robust')
xlabel('FR to punish cue');
ylabel('RT to punish cue');
fnm = fullfile(resdir,'FR_vs_RT_2.fig'); % save figures
saveas(gcf,fnm);
fnm = fullfile(resdir,'FR_vs_RT_2.jpg');
saveas(gcf,fnm);
set(gcf, 'renderer', 'painters')
fnm = fullfile(resdir,'FR_vs_RT_2.eps');
saveas(gcf,fnm);

% Regression of lick rate and firing rate
polypredcicall(FR1_avg,LR1_avg,0.95,'robust')
xlabel('FR to reward cue');
ylabel('LR to reward cue');
fnm = fullfile(resdir,'FR_vs_LR_1.fig'); % save figures
saveas(gcf,fnm);
fnm = fullfile(resdir,'FR_vs_LR_1.jpg');
saveas(gcf,fnm);
set(gcf, 'renderer', 'painters')
fnm = fullfile(resdir,'FR_vs_LR_1.eps');
saveas(gcf,fnm);

polypredcicall(FR2_avg,LR2_avg,0.95,'robust')
xlabel('FR to punish cue');
ylabel('LR to punish cue');
fnm = fullfile(resdir,'FR_vs_LR_2.fig'); % save figures
saveas(gcf,fnm);
fnm = fullfile(resdir,'FR_vs_LR_2.jpg');
saveas(gcf,fnm);
set(gcf, 'renderer', 'painters')
fnm = fullfile(resdir,'FR_vs_LR_2.eps');
saveas(gcf,fnm);

function avg_LR = calc_avg_LR(Licks, StimOn)
ntrials = length(StimOn);
LR = nan(1,ntrials);
for n = 1:ntrials
    lick_ts = Licks{n};
    LR(n) = length(lick_ts((lick_ts > StimOn(n)) & (lick_ts < StimOn(n)+1)));
end
avg_LR = nanmean(LR);