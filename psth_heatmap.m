function psth_heatmap(psth1,psth2,time,resdir)
%PSTH_HEATMAP   Heatmap for all PSTHs.
%   PSTH_HEATMAP(PSTH1,PSTH2,TIME,RESDIR) plots sorted PSTHs color coded.
%   Plots are saved in RESDIR.
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

% Time windows
twin = [-0.2 0.5];  % for sorting
sr = 1000;
tinx = time >= twin(1) * sr & time <= twin(2) * sr;

plotwin = [-3 3];  % for plotting
pwinx = time >= plotwin(1) * sr & time <= plotwin(2) * sr;

% Sort PSTHs (based on max value in the second one)
mv1 = max(psth2(:,tinx),[],2);
[~, ia] = sort(mv1,'descend');

% Plot heatmaps
NumCells = size(psth1,1);
clims = prctile(psth2(:),[20 99.5]);  % color range

H1 = figure;
imagesc(time(pwinx),1:NumCells,psth1(ia,pwinx),clims);  % first
colormap hot
axis square

H2 = figure;
imagesc(time(pwinx),1:NumCells,psth2(ia,pwinx),clims);   % second
colormap hot
axis square

% Save
fnm = fullfile(resdir,'psth_heatmap_T1.fig');
fnm2 = fullfile(resdir,'psth_heatmap_T2.fig');
saveas(H1,fnm);  % save first heatmap
close(H1)
saveas(H2,fnm2);  % save second heatmap
close(H2)