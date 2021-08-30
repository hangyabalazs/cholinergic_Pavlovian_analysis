function psth_by_label(psth,time,labels,resdir)
%PSTH_BY_LABEL   Group PSTHs by label.
%   PSTH_BY_LABEL(PSTH,TIME,LABELS,RESDIR) groups PSTHs by LABELS and plots
%   them as a function of TIME, both individually and averaged. The
%   resulting plots are saved in RESDIR.
%
%   See also AVG_PSTH_CHOLINERGIC.

%   Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.hu
%   12-21-2020

% Results directory
if ~isfolder(resdir)
    mkdir(resdir)
end

% Grouping
ulabels = unique(labels);  % unique labels
NumLabels = length(ulabels);  % number of labels
inx_psth = cell(1,NumLabels);  % indices of PSTHs per label
NumPSTH = nan(1,NumLabels);  % number of PSTHs per label
for iL = 1:NumLabels
    inx_psth{iL} = find(strcmp(ulabels{iL},labels));
    NumPSTH(iL) = length(inx_psth{iL});
end
[~, ia] = sort(NumPSTH,'descend');
ulabels = ulabels(ia);  % sort in descending order of number of PSTHs for better plotting
NumPSTH = NumPSTH(ia);
inx_psth = inx_psth(ia);

% Plot
H1 = figure;  % all PSTHs
hold on
H2 = figure;  % average PSTH per label
hold on
[lss1, lss2] = deal(nan(1,NumLabels));
for iL = 1:NumLabels
    lpsth = psth(inx_psth{iL},:);
    for iN = 1:NumPSTH(iL)
        c = [iL/NumLabels 1-iL/NumLabels iL/NumLabels];
        figure(H1)
        L = plot(time,lpsth(iN,:),'Color',c,'LineWidth',2);
    end
    lss1(iL) = L;
    figure(H2)
    L = plot(time,mean(lpsth,1),'Color',c,'LineWidth',2);
    lss2(iL) = L;
end
figure(H1)
legend(lss1,ulabels)
figure(H2)
legend(lss2,ulabels)

% Save
fnm = fullfile(resdir, 'all_spth.fig');
fnm2 = fullfile(resdir, 'average_psth.fig');
saveas(H1,fnm);  % save all PSTHs plot
close(H1)
saveas(H2,fnm2);  % save average PSTHs plot
close(H2)