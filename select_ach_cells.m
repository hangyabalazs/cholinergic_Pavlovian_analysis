function achcells = select_ach_cells(achcells)
%SELECT_ACH_CELLS   Select cholinergic neurons from CellBase.
%   SELECT_ACH_CELLS(ACHCELLS) selects cells that are assigned to
%   cholinergic ('ChAT+' property set to 1) from CELLIDS (default, all cell
%   IDs in CellBase). Only well-isolated units (ID > 20, L_ratio < 0.15)
%   are returned.
%
%   See also CHOLINERGIC_ANALYSIS_PAVLOVIAN, VPSELECTCELLS.

%   Balazs Hangya and Panna Hegedus
%   Laboratory of Systems Neuroscience
%   Institute of Experimental Medicine, Budapest, Hungary

%   Code review: 5/9/20

% Input arguments
narginchk(0,1)
if nargin < 1
    achcells = [];
end

% Load CellBase
load(getpref('cellbase','fname'));

% List of cellIDs
if isempty(achcells)
    Lratio = getvalue('Lr_PC');
    ID = getvalue('ID_PC');
    selstr = '"ChAT+"==1';    % select cholinergic cells
    achcells = selectcell(selstr);
else
    if isnumeric(achcells)  % 'vpcells' can be index set or list of cell IDs
        achcells = CELLIDLIST(achcells);
    else
        if ischar(achcells)
            achcells = {achcells};   % only one cell ID
        end
    end
end
achcells = achcells(:)';   % convert to row vector