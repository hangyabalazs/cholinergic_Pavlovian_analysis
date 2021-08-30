function ach_acg(cellids, resdir, issave)
%ACH_ACG   Auto-correlation analysis.
%   ACH_ACG calculates auto-correlations for cholinergic neurons in CELLIDS in 500 ms windows.
%   ACG results are saved to RESDIR if ISSAVE is true.
%
%   See also ACG.

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   16-04-2020

% Code review: BH 5/11/20

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
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

% ACG
segfilter = 'stim_excl_vp'; %should I write a different filter? This was for VP but its the same for the cholinergic neurons as well
filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
acg(achcells,0.5,'resdir',resdir,...
     'segfilter',segfilter,'filterinput',filterinput,...
     'minspikeno',100,'maxspikeno',10000,'issave',issave);