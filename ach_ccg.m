function ach_ccg(cellids,resdir,issave)
%ACH_CCG   Cross-correlation analysis.
%   ACH_CCG calculates cross-correlations for BF cholinergic neurons listed
%   in CELLIDS in 50 ms windows. CCG results are saved for non-tetrode 
%   pairs and tetrode pairs separately to RESDIR if ISSAVE is true.
%
%   See also CCG.

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   16-Apr-2020

%   Code review: BH 5/11/20

% Directories
resdir_tetrode = fullfile(resdir,'tetrodepairs');   % results directory for tetrode pairs
resdir_nontetrode = fullfile(resdir,'nontetrodepairs');   % results directory for non-tetrode pairs
if ~isfolder(resdir_tetrode)
    mkdir(resdir_tetrode)
end
if ~isfolder(resdir_nontetrode)
    mkdir(resdir_nontetrode)
end

% Input argument check
narginchk(0,3);
if nargin < 3
    issave = true;   % default saving behavior
end
if nargin < 1
    achcells = select_ach_cells([]);
else
    achcells = cellids;
end
achcells = achcells';

% CCG
segfilter = 'stim_excl_vp';
filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
ccg(achcells,0.05,'whichcells','nontetrodepairs','resdir',resdir_nontetrode,...
     'segfilter',segfilter,'filterinput',filterinput,...
     'minspikeno',100,'maxspikeno',10000,'issave',issave);   % non-tetrode pairs
ccg(achcells,0.05,'whichcells','tetrodepairs','resdir',resdir_tetrode,...
    'segfilter',segfilter,'filterinput',filterinput,...
    'minspikeno',100,'maxspikeno',10000,'issave',issave);   % tetrode pairs