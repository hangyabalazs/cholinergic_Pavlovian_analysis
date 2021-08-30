function [g1, g2] = countlicks(animalID, sessionID)
%COUNTLICKS Behavioral discrimination of cues in Pavlovian conditioning.
%   [LR1, LR2] = COUNTLICKS(ANIMALID,SESSIONID) calculates number of licks
%   in a 1.2 seconds window after stimulus onset for two trial types in
%   Pavlovian conditioning. Lick counts are calculated for the session
%   specified by ANIMALID and SESSIONID and returned in LR1 and LR2 for the
%   two trial types.
%
%   To calculate lick rates (lick counts normalized by elapsed time), use
%   this example:
%   [TrialType1, TrialType2] = countlicks(animalID, sessionID);   % licks per trial type for each trial
%   T1 = nansum(TrialType1)/(sum(~isnan(TrialType1))*1.2);   
%   T2 = nansum(TrialType2)/(sum(~isnan(TrialType2))*1.2);
%
%   See also ANTICIPATORY_STATS.

%   Panna Hegedus, 05-Feb-2020
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu 

% Code review: BH 4/8/20

% Load behavior data
cbdir = getpref('cellbase','datapath');
datapath = fullfile(cbdir,animalID,sessionID,'TrialEvents.mat');
VE = load(datapath);
if length(unique(VE.TrialType)) == 2
    ntrials = length(VE.NTrials); % length of each group
    
    % the frist 20 trials can be free access of water in the
    % training protocol, therefore these trials are skipped
    [g1, g2] = deal(zeros(1,ntrials));   % preallocate space
    for cT = 21:ntrials  % loop through trials from (skip first 20 trials)
        ws = VE.StimulusOn(cT);
        we = ws + 1.2; % window: stimulus duration + 200ms delay
        if VE.TrialType(cT) ==1 && (~isempty(VE.LickIn{cT})) % for TrialType = 1
            licks = VE.LickIn{cT}; % lick
            licks = licks(licks>ws);
            licks = licks(licks<we);
            g1(1, cT) = length(licks);
        elseif VE.TrialType(cT) ==2 && (~isempty(VE.LickIn(cT))) % for TrialType = 2
            licks = VE.LickIn{cT};
            licks = licks(licks>ws);
            licks = licks(licks<we);
            g2(1, cT) = length(licks);
        end
    end
else
    g1 = NaN;
    g2 = NaN; % do not make statistics, if there is only 1 trialtype
end
end