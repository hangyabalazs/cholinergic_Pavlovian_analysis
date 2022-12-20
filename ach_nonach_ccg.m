function ach_nonach_ccg(cellids,resdir,issave)
%ACH_NONACH_CCG(CELLIDS, RESDIR, ISSAVE)   Cross-correlation analysis of cholinergic and noncholinergic neurons.
%   ACH_NONACH_CCG calculates cross-correlations for BF cholinergic neurons
%   (listed in CELLIDS) and noncholinergic neurons recorded in the same
%   session in 50 ms windows. CCG results are saved to RESDIR if ISSAVE is true.
%
%   See also CCG.

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   19-May-2020

resdir_allpairs = fullfile(resdir,'ach_nonach_pairs');   % results directory for tetrode pairs

if ~isfolder(resdir_allpairs) % make results directory if non existent
    mkdir(resdir_allpairs)
end

% Input argument check
narginchk(0,3);
if nargin < 3
    issave = true;   % default saving behavior
end
if nargin < 1 % if no input - cholinergic cells are chosen automatically
    achcells = select_ach_cells([]);
else
    achcells = cellids; %cholinergic cells
end
achcells = achcells';

% CCG
segfilter = 'stim_excl_vp'; %filter for ccg
filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};

aID = unique(getvalue('RatID_tag', achcells)); % animal IDs

for i = 1:length(aID)
    aID_achcell = findachcell(aID{i});    % find cholinergic cell of given animal ID
    sID = getvalue('SessionID_tag', aID_achcell); % find sessionID of the cell
    for n = 1:length(sID)
    session_cellids = findcell('rat',aID{i}, 'session', sID{n}); % find all other cells in the session where the cholinergic neuron was recorded
    nonach_cellids = setdiff(session_cellids, achcells); % separate the noncholinergic ones
    session_achcells = intersect(achcells, session_cellids); % separate the cholinergic ones - there can be multiple cells
    
    for k = 1:length(session_achcells)
        cAchcell = session_achcells{k}; % choose a cholinergic cell
        for j = 1:length(nonach_cellids)
            cellpairs = {cAchcell nonach_cellids{j}}; % make a cell pair from a cholinergic and a noncholinergic cell
            ccg(cellpairs,0.05,'whichcells','allpairs','resdir',resdir_allpairs,...
                'segfilter',segfilter,'filterinput',filterinput,...
                'minspikeno',100,'maxspikeno',10000,'issave',issave);   % cholinergic-non cholinergic pairs
        end
    end
end
end

% -------------------------------------------------------------------------
function aID_achcell = findachcell(aID)
% Load CellBase
load(getpref('cellbase','fname')); %tried to do it with 'selectcell' bbut could not succeed
ischat = getvalue('ChAT+');
tag = getvalue('RatID_tag');
selstr = ischat==1 & strcmp(tag, aID);    % select cholinergic cells
aID_achcell = CELLIDLIST(selstr);





