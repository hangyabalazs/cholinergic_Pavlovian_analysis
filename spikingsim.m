function spikingsim(cellids,model_values)
%SPIKINGSIM   Simulate spiking bassed on TDRL model.
%   SPIKINGSIM(CELLIDS,ETA1,ETA2,MODEL_VALUES) generates a virtual,
%   simulated spike train for cells in CELLIDS. A frequency-matched Poisson
%   process is summed with evoked spikes, where the number of spikes is
%   based on simulated values of a best-fit TDRL model and simulated spike
%   times are drawn from a Gaussian distribution with fixed latency and 
%   jitter. The virtual spike train is saved in the respective CellBase
%   folder with a VV tag, using VV_<tetrode>_<unit> convention.
%
%   See also CHOLINERGIC_RFMODEL_MAIN.

%   Balazs Hangya, 20-Nov-2020
%   Institute of Experimental Medicine
%   hangya.balazs@koki.hu 

% Loop through cells
NumCells = length(cellids);
for iC = 1:NumCells
    cellid = cellids{iC};
    model_value = model_values{iC};
    main(cellid,model_value)
end

% -------------------------------------------------------------------------
function main(cellid,model_value)

% Load spike times
SP = loadcb(cellid,'Spikes');
datalength = SP(end) - SP(1);
FR = (length(SP) - 1) / datalength;

% Load trial events
VE = loadcb(cellid,'TrialEvents');   % load events

% Parameters
mx = datalength * 1000;   % data length in ms
baselinerate = FR;   % baseline spike rate (Poisson) in Hz
win_cue = 0.5;  % window size for cue (s)
win_reinf = 0.2;  % window size for reinforcement (s)

% Response latencies
esp_lat_cue = 0.105; % latency of evoked spikes after cue in s (Gaussian)
esp_jit_cue = 0.001; % jitter of evoked spikes after cue in s (Gaussian)
esp_lat_rew = 0.092; % latency of evoked spikes after reward in s (Gaussian)
esp_jit_rew = 0.001; % jitter of evoked spikes after reward in s (Gaussian)
esp_lat_pun = 0.075; % latency of evoked spikes after punishment in s (Gaussian)
esp_jit_pun = 0.001; % jitter of evoked spikes after punishment in s (Gaussian)

% Number of evoked spikes based on TDRL model
fresp_T1_cue = model_value(1,1) - FR;     % eveoked spike rate for cues in type1 trials
fresp_T2_cue = model_value(1,2) - FR;     % eveoked spike rate for cues in type2 trials
fresp_T1_rew = model_value(2,1) - FR;     % eveoked spike rate for reward in type1 trials
fresp_T2_rew = model_value(2,2) - FR;     % eveoked spike rate for reward in type2 trials
fresp_T1_pun = model_value(3,1) - FR;     % eveoked spike rate for punishment in type1 trials
fresp_T2_pun = model_value(3,2) - FR;     % eveoked spike rate for punishment in type2 trials
nesp_T1_cue = round(fresp_T1_cue / win_cue);     % number of eveoked spikes for cues in type1 trials
nesp_T2_cue = round(fresp_T2_cue / win_cue);     % number of eveoked spikes for cues in type2 trials
nesp_T1_rew = round(fresp_T1_rew / win_reinf);     % number of eveoked spikes for reward in type1 trials
nesp_T2_rew = round(fresp_T2_rew / win_reinf);     % number of eveoked spikes for reward in type2 trials
nesp_T1_pun = round(fresp_T1_pun / win_reinf);     % number of eveoked spikes for punishment in type1 trials
nesp_T2_pun = round(fresp_T2_pun / win_reinf);     % number of eveoked spikes for punishment in type2 trials

% Simulate event train 
type1inx = VE.TrialType == 1;  % indices for Type1 trials
E_T1_cue = VE.TrialStart(type1inx) + VE.StimulusOn(type1inx);   % cue presentation times in Type1 trials
NumT1Cue = length(E_T1_cue); % number of Type1 trials
type2inx = VE.TrialType == 2;  % indices for Type2 trials
E_T2_cue = VE.TrialStart(type2inx) + VE.StimulusOn(type2inx);   % cue presentation times in Type2 trials
NumT2Cue = length(E_T2_cue); % number of Type2 trials
rewinxT1 = VE.AllReward == 1; % indices for rewarded type1 trials
E_T1_rew = VE.TrialStart(rewinxT1) + VE.DeliverFeedback(rewinxT1);   % reward delivery times in Type1 trials
NumT1Rew = length(E_T1_rew); % number of rewards in Type1 trials
rewinxT2 = VE.AllReward == 2; % indices for rewarded Type2 trials
E_T2_rew = VE.TrialStart(rewinxT2) + VE.DeliverFeedback(rewinxT2);   % reward delivery times in Type2 trials
NumT2Rew = length(E_T2_rew); % number of rewards in Type2 trials
puninxT1 = VE.Punishment == 1; % indices for punished Type1 trials
E_T1_pun = VE.TrialStart(puninxT1) + VE.DeliverFeedback(puninxT1);   % punishment delivery times in Type1 trials
NumT1Pun = length(E_T1_pun); % number of punishments in Type1 trials
puninxT2 = VE.Punishment == 2; % indices for punished Type2 trials
E_T2_pun = VE.TrialStart(puninxT2) + VE.DeliverFeedback(puninxT2);   % punishment delivery times in Type2 trials
NumT2Pun = length(E_T2_pun); % number of punishments in Type2 trials

% Simulate spike train
S0 = randpoisson(mx/1000*baselinerate,mx)' / 1000;    % generate background Poisson spiking
S1 = repmat(E_T1_cue,nesp_T1_cue,1) + esp_lat_cue + randn(nesp_T1_cue,NumT1Cue) * esp_jit_cue;  % generate 'cue-evoked' spikes in Type1 trials
% S1 = repmat(E_T1_cue,nesp_T1_cue*100,1)+2;
S1 = sort(S1(:));
S2 = repmat(E_T2_cue,nesp_T2_cue,1) + esp_lat_cue + randn(nesp_T2_cue,NumT2Cue) * esp_jit_cue;  % generate 'cue-evoked' spikes in Type2 trials
S2 = sort(S2(:));
S3 = repmat(E_T1_rew,nesp_T1_rew,1) + esp_lat_rew + randn(nesp_T1_rew,NumT1Rew) * esp_jit_rew;  % generate 'reward-evoked' spikes in Type1 trials
S3 = sort(S3(:));
S4 = repmat(E_T2_rew,nesp_T2_rew,1) + esp_lat_rew + randn(nesp_T2_rew,NumT2Rew) * esp_jit_rew;  % generate 'reward-evoked' spikes in Type2 trials
S4 = sort(S4(:));
S5 = repmat(E_T1_pun,nesp_T1_pun,1) + esp_lat_pun + randn(nesp_T1_pun,NumT1Pun) * esp_jit_pun;  % generate 'punishment-evoked' spikes in Type1 trials
S5 = sort(S5(:));
S6 = repmat(E_T2_pun,nesp_T2_pun,1) + esp_lat_pun + randn(nesp_T2_pun,NumT2Pun) * esp_jit_pun;  % generate 'punishment-evoked' spikes in Type2 trials
S6 = sort(S6(:));
VS = sort([S0; S1; S2; S3; S4; S5; S6]);   % all spike times, in s

% Save simulated spike train
[animalID, sessionID, Tetrode, Unit] = cellid2tags(cellid);
resdir = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
fnm = fullfile(resdir,['VV' num2str(Tetrode) '_' num2str(Unit)]);
save(fnm,'VS')

% Pre-align spikes
prealignVSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_pavlovian,'filetype','event','ifsave',1,'ifappend',0)

% Diagnostic CCG
ncc1 = E_T1_cue(:);
ncc2 = VS(:);
wn = 2.7;
% diagnosticCCG(ncc1,ncc2,wn)

% -------------------------------------------------------------------------
function diagnosticCCG(ncc1,ncc2,wn)

% Calculate spike times in milliseconds
sr = 1000;
nc1 = ncc1 * sr;
nc2 = ncc2 * sr;
mn = min(nc1(1),nc2(1));  % only relative spike times count; avoid out of memory
nc1 = nc1 - mn;
nc2 = nc2 - mn;
nc1(nc1<0.5) = [];  % drop spikes before 0.5 ms(!) to avoid error in line 39
nc2(nc2<0.5) = [];
wn2 = wn * sr;  % window size in ms

% Crosscorrelogram
zunit1 = zeros(1,round(max([nc1; nc2]))+5);
zunit2 = zunit1;
zunit1(round(nc1)) = 1;
zunit2(round(nc2)) = 1;
[ccr, lags] = xcorr(zunit2,zunit1,wn2);     % 1->2; window: -wn ms - wn ms

% Plot
H1 = figure;
bar(lags,ccr,'FaceColor','black')