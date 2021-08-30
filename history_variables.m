function history_variables(cellids)
%HISTORY_VARIABLES   Save outcome history variables in TrialEvents.
%   HISTORY_VARIABLES(CELLIDS) saves a 'PreviousOutcome' variable in
%   TrialEvents that indicates if previous trials was rewarded (1),
%   punished (2) or ended with outcome omission (3).
%
%   See also AVG_PSTH_CHOLINERGIC.

%   Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu

% Create history variables for all cells
NumCells = length(cellids);
for iC = 1:NumCells
    
    % Previous outcome variables
    TE = loadcb(cellids{iC},'TrialEvents');
    NumTrials = length(TE.NTrials); % number of trials
    if ~isfield(TE,'PreviousOutcome')  % add lick variables to TrialEvents
        previous_outcome_varaibles(cellids{iC},TE,NumTrials)
    end
end

% -------------------------------------------------------------------------
function previous_outcome_varaibles(cellid,TE,NumTrials)

TE.PreviousOutcome = nan(1,NumTrials);
for iT = 2:NumTrials
    if ~isnan(TE.Omission(iT-1))  % check previous trial
        TE.PreviousOutcome(iT) = 3;  % omission
    elseif ~isnan(TE.RewardedTrials(iT-1))
        TE.PreviousOutcome(iT) = 1;  % reward
    elseif ~isnan(TE.PunishedTrials(iT-1))
        TE.PreviousOutcome(iT) = 2;  % punishment
    end
end
aID = getvalue('RatID_tag', cellid);
sID = getvalue('SessionID_tag', cellid);
savedir = fullfile(getpref('cellbase', 'datapath'), aID{1}, sID{1});
savename = 'TrialEvents.mat'; % save TE struct
savename2 = 'TE.mat'; % save TE struct
save(fullfile(savedir, savename),'-struct','TE');
save(fullfile(savedir, savename2),'-struct','TE');