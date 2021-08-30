function [a1, a2, S, PARAM, fr, err] = best_fit_cholinergic(cellid,spike_data,eps,varargin)
%BEST_FIT_CHOLINERGIC   TDRL model fit on cholinergic spiking data.
%   [A, S] = BEST_FIT_CHOLINERGIC(CELLID, SPIKE_DATA, EPS) generates a
%   contingency matrix of based on trial parameters used for model fitting.
%   Each row corresponds to a cue or a reinforcement event (concatenated)
%   and columns contain reward expectation, punishment expectation,
%   indicator variables for reward and punishment, and SD of spiking data
%   for the noise model. Finally, ACH_MODEL is called for model fitting.
%
%   See also CHOLINERGIC_RFMODEL_MAIN and ACH_MODEL.

%   Panna Hegedus
%   Insititute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   15-Apr-2020

%   Code review: BH 4/16/20

% Parsing inputs
prs = inputParser;
addRequired(prs,'cellid',@(s)isempty(s)|iscell(s)|iscellstr(s)|ischar(s)) % cellIDs - cells
addRequired(prs,'spike_data',@iscell)  % firing rate data
addRequired(prs,'eps',@isnumeric)   % SD of firing rate data for the moise model
addParameter(prs,'isreversed',false,@logical)   % if true, run with reverse contingencies as control
parse(prs,cellid,spike_data,eps,varargin{:})
g = prs.Results;

% Trial events
VE = loadcb(cellid,'TrialEvents');   % load events

% Parameters
numTrials = VE.NTrials(1);
E_R = {0.8, 0.25};  % task contingencies - reward probability in TrialType 1 and 2
E_P = {0.1, 0.65};  % task contingencies - punishment probability in TrialType 1 and 2
if g.isreversed
    E_R = {0.1, 0.65};   % run the model with reverse contingencies, as control
    E_P = {0.8, 0.25};
end
P = {0, 1};
R = {0, 1};

% Find omission trials
ominx = (VE.Omission==1 | VE.Omission==2);

% Reward expectation
ER_cue = nan(1,numTrials);   % reward expectation after the cue
ER_reinf = nan(1,numTrials);   % reward expectation at reinforcement
ER_cue(VE.TrialType==1) = E_R{1};
ER_cue(VE.TrialType==2) = E_R{2};
ER_reinf(VE.AllReward==1 | VE.Punishment==1) = E_R{1};  % rewarded or punished in TrialType1 (likely reward)
ER_reinf(VE.AllReward==2 | VE.Punishment==2) = E_R{2};  % rewarded or punished in TrialType2 (likely punishment)
ER = [ER_cue(~ominx) ER_reinf(~ominx)]';  % concatenate cues and reinforcements in all trials except omissions

% Punishment expectation
EP_cue = nan(1,numTrials);   % punishment expectation after the cue
EP_reinf = nan(1,numTrials);   % punishment expectation at reinforcement
EP_cue(VE.TrialType==1) = E_P{1};
EP_cue (VE.TrialType==2) = E_P{2};
EP_reinf(VE.AllReward==1 | VE.Punishment==1) = E_P{1};  % rewarded or punished in TrialType1 (likely reward)
EP_reinf(VE.AllReward==2 | VE.Punishment==2) = E_P{2};  % rewarded or punished in TrialType2 (likely punishment)
EP = [EP_cue(~ominx) EP_reinf(~ominx)]';  % concatenate cues and reinforcements in all trials except omissions

% Reward
R_reinf = nan(1,numTrials);  % indicator for reward
R_cue = zeros(1,numTrials);
R_reinf(VE.AllReward==1 | VE.AllReward==2) = R{2}; 
R_reinf(VE.Punishment==1 | VE.Punishment==2) = R{1};
R_all = [R_cue(~ominx) R_reinf(~ominx)]';

% Punishment
P_reinf = nan(1,numTrials);   % indicator for punishment
P_cue = zeros(1,numTrials);
P_reinf(VE.AllReward==1 | VE.AllReward==2) = P{1}; 
P_reinf(VE.Punishment==1 | VE.Punishment==2) = P{2};
P_all = [P_cue(~ominx) P_reinf(~ominx)]';

% Epsilon
E_reinf = nan(1,numTrials);
E_cue = repmat(eps(1),1,numTrials);   % SD of spike no/firing rate after the cue
E_reinf(VE.AllReward==1 | VE.AllReward==2) = eps(2);    % SD of spike no/firing rate after reward
E_reinf(VE.Punishment==1 | VE.Punishment==2) = eps(3);   % SD of spike no/firing rate after punishment
E = [E_cue(~ominx) E_reinf(~ominx)]';

% All trial information
PARAM = [ER EP R_all P_all E];

% Spikes
fr2 = nan(numTrials,1);
rinx = VE.AllReward==1 | VE.AllReward==2;   % rewarded trials
pinx = VE.Punishment==1 | VE.Punishment==2;   % punished trials

fr1 = spike_data{1}';   % spike no/FR after cue
fr2(rinx,1) = spike_data{2}';   % spike no/FR after reward
fr2(pinx,1) = spike_data{3}';   % spike no/FR after punishment
fr1 = fr1(~ominx);
fr2 = fr2(~ominx);
fr = [fr1; fr2];    % concatenate cue and reinforcement 

% Fit TDRL model
[a1, a2, S, err] = Ach_model(PARAM, fr);