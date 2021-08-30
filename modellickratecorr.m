function modellickratecorr(cellids,eta1,eta2,resdir)
%MODELLICKRATECORR   Correlation between model parameters and behavioral discrimination.
%   MODELLICKRATECORR(CELLIDS,ETA1,ETA2) calculates correlation between
%   lick rate difference between trial types in a given session (see
%   COUNTLICKS) and TDRL model parameters ETA1 and ETA2 (see
%   CHOLINERGIC_RFMODEL_MAIN).
%
%   See also COUNTLICKS and CHOLINERGIC_RFMODEL_MAIN.

%   Balazs Hangya, 20-Nov-2020
%   Institute of Experimental Medicine
%   hangya.balazs@koki.hu 

% Results directory
if ~isfolder(resdir)
    mkdir(resdir)
end

% Animal and session IDs
animals = getvalue('RatID_tag',cellids);
sessions = getvalue('SessionID_tag',cellids);

% Calculate lick rates
NumCells = length(sessions);
[T1, T2] = deal(nan(NumCells,1));
for iC = 1:NumCells  % loop through cells
    animalID = animals{iC};
    sessionID = sessions{iC};
    [TrialType1, TrialType2] = countlicks(animalID, sessionID);   % lick counts per trial type for each trial
    T1(iC) = nansum(TrialType1)/(sum(~isnan(TrialType1))*1.2);   % lick rate, TrialType1
    T2(iC) = nansum(TrialType2)/(sum(~isnan(TrialType2))*1.2);   % lick rate, TrialType2
end

% Calculate regression
polypredcicall(T1-T2,eta1,0.95,'robust')
xlabel('anticipatory lick rate difference');
ylabel('eta_1');
fnm = fullfile(resdir,'lickratediff_vs_eta1.fig'); % save figures
saveas(gcf,fnm);
fnm = fullfile(resdir,'lickratediff_vs_eta1.jpg');
saveas(gcf,fnm);
set(gcf, 'renderer', 'painters')
fnm = fullfile(resdir,'lickratediff_vs_eta1.eps');
saveas(gcf,fnm);

polypredcicall(T1-T2,eta2,0.95,'robust')
xlabel('anticipatory lick rate difference');
ylabel('eta_2');
fnm = fullfile(resdir,'lickratediff_vs_eta2.fig'); % save figures
saveas(gcf,fnm);
fnm = fullfile(resdir,'lickratediff_vs_eta2.jpg');
saveas(gcf,fnm);
set(gcf, 'renderer', 'painters')
fnm = fullfile(resdir,'lickratediff_vs_eta2.eps');
saveas(gcf,fnm);

polypredcicall(T1-T2,eta1+eta2,0.95,'robust')
xlabel('anticipatory lick rate difference');
ylabel('eta_1 + eta_2');
fnm = fullfile(resdir,'lickratediff_vs_etasum.fig'); % save figures
saveas(gcf,fnm);
fnm = fullfile(resdir,'lickratediff_vs_etasum.jpg'); % save figures
saveas(gcf,fnm);
set(gcf, 'renderer', 'painters')
fnm = fullfile(resdir,'lickratediff_vs_etasum.eps'); % save figures
saveas(gcf,fnm);
end