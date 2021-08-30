function cholinergic_PSTH(area, resdir)
%   CHOLINERGIC_PSTH(AREA, RESDIR) non-adaptive, adaptive and doubly
%   adaptive PSTH of cholinergic BF neurons
%   CHOLINERGIC_PSTH generates non-adaptive, adaptive and doubly adaptive
%   PSTH of BF cholinergic neurons of specified location (AREA). Spikes are aligned to stimulus onset or
%   feedback delivery. Neural activity is filtered to different cues and
%   outcomes (reward and punishment) as well. Results are saved to RESDIR.
%
%   See also VIEWCELL2B.
%
%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   18-May-2020



if ~isfolder(resdir) % make results directory if it does not exists
    mkdir(resdir);
end

narginchk(0,2)
if nargin < 2
    resdir = 'F:\pavlovian_cholinergic_cellbase\_adaptive_PSTH_all';
end


selstr = [' "ChAT+"==1&strcmp("Area1", {[''' area ''']})'];
cellids = selectcell(selstr);   % select cholinergic neurons

tags = {'nonadaptive' 'adaptive' 'doubly_adaptive'} % isadaptive
alignment = {{'DeliverAllFeedback' 'Feedback'}}; %triggername and filter pairs for viewcell2b

for i = 1:length(tags) %loop through tags
    
    tag = tags{i};
    switch tag
        case 'nonadaptive'
            a = 0;
        case 'adaptive'
            a = 1;
        case 'doubly_adaptive'
            a = 2;
    end
    for j = 1:length(alignment)
        trigger = alignment{j}{1};
        filter = alignment{j}{2};
        for k = 1:length(cellids)
            pause(0.01)
            viewcell2b(cellids{k},'TriggerName',trigger,'SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions',['#' filter],'isadaptive', a, 'window',[-3 3])
            
            cellidt = cellids{k};
            cellidt(cellidt=='.') = '_';
            
            fnm = fullfile(resdir,[cellidt '_' area '_' filter '_' tag '.jpeg']);   % save
            fnm = fullfile(resdir,[cellidt  '_' filter '_' tag '.jpeg']);
            saveas(gcf,fnm)
            close(gcf)
        end
    end
end