function cue_reward_correlation_cholinergic(cellids, resdir)
%CUE_REWARD_CORRELATION_CHOLINERGIC   Correlation between reward and cue
%response of cholinergic neurons.
%   CUE_REWARD_CORRELATION_CHOLINERGIC(CELLIDS, RESDIR) plots a scatter
%   plot and performs linear regression on cue and reward response of
%   cholinergic neurons. Correlation coefficient is calculated as well.
%   Correlations are examined within each cell and across animals as well.
%   Results are saved to RESDIR.
%
%   See also ULTIMATE_PSTH, REGRESSIONPLOT, SCATTER and CORRCOEF.

%   Panna Hegedus
%   Institute of Experimental Medicine
%   panna.hegedus@koki.mta.hu
%   16-Nov-2020

%   Code review: BH 12/2/20

% Results directory
if ~isfolder(resdir)
    mkdir(resdir)
end

% Cells
if nargin <1
    cellids = select_ach_cells([]);
end

% Time window
wn = [-4 4];
dt = 0.001;
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector

twin1 = [0 0.2]; % test window for reward response
twin2 = [0 0.5]; % test window for cue response

bwin1 = [-3 -2]; % baseline window for reward response
bwin2 = [-1 0]; % baseline window for cue response

% Index for time 0
st = abs(wn(1)) / dt;   % in ms
nullindex = st + 1;

% PSTH
R1 = runanalysis(@ultimate_psth,...   % Reward response for trial type 1
    'trial', 'DeliverAllFeedback', wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
    'event_filter','custom', 'filterinput','AllReward==1','maxtrialno',Inf,...
    'baselinewin',bwin1,'testwin',twin1,'relative_threshold',0.1,'cellids',cellids);
R2 = runanalysis(@ultimate_psth,...   % Reward response for trial type 2
    'trial', 'DeliverAllFeedback', wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
    'event_filter','custom', 'filterinput','AllReward==2','maxtrialno',Inf,...
    'baselinewin',bwin1,'testwin',twin1,'relative_threshold',0.1,'cellids',cellids);

C1 = runanalysis(@ultimate_psth,...   % Cue response for trial type 1
    'trial', 'StimulusOn', wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
    'event_filter','custom', 'filterinput','AllReward==1','maxtrialno',Inf,...
    'baselinewin',bwin2,'testwin',twin2,'relative_threshold',0.1,'cellids',cellids);
C2 = runanalysis(@ultimate_psth,...   % Cue response for trial type 2
    'trial', 'StimulusOn', wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
    'event_filter','custom', 'filterinput','AllReward==2','maxtrialno',Inf,...
    'baselinewin',bwin2,'testwin',twin2,'relative_threshold',0.1,'cellids',cellids);

A = cellfun(@cellid2tags, cellids, 'UniformOutput', false); % get animal IDs
animals = unique(A);
NumAnimals = length(animals);

WNt_reward = [twin1(1)/dt+nullindex twin1(2)/dt+nullindex-1];   % test window; convert to indices
WNt_cue = [twin2(1)/dt+nullindex twin2(2)/dt+nullindex-1];   % test window; convert to indices

WNb_reward = [bwin1(1)/dt+nullindex bwin1(2)/dt+nullindex-1]; % baseline window; convert to indices
WNb_cue = [bwin2(1)/dt+nullindex bwin2(2)/dt+nullindex-1];

BR_R1 = R1(:,5); % binrasters
BR_R2 = R2(:,5);
BR_C1 = C1(:,5);
BR_C2 = C2(:,5);

% Pool data across cells within mice
[cue_per_animal1,cue_per_animal2,reward_per_animal1, reward_per_animal2]  = deal(nan(1,NumAnimals)); % preallocate
for k = 1:NumAnimals % loop through animals
    inx = ~cellfun(@isempty, (strfind(cellids, animals{k}))); % indices of cells for the current animal
    
    r1 = vertcat(BR_R1{inx}); % concatenate reward responses
    r2 = vertcat(BR_R2{inx});
    
    c1 = vertcat(BR_C1{inx}); % concatenate cue responses
    c2 = vertcat(BR_C2{inx});
    
    [cue_per_cell1, cue_per_cell2, reward_per_cell1, reward_per_cell2] = deal([]);
    
    % Get trial-by-trial cue and reward responses (baseline substracted)
    c1_ind = sum(c1(:,WNt_cue(1):WNt_cue(2)),2)/twin2(2)-sum(c1(:,WNb_cue(1):WNb_cue(2)),2)/abs(bwin2(2)-bwin2(1));
    c2_ind = sum(c2(:,WNt_cue(1):WNt_cue(2)),2)/twin2(2)-sum(c2(:,WNb_cue(1):WNb_cue(2)),2)/abs(bwin2(2)-bwin2(1));
    
    r1_ind = sum(r1(:,WNt_reward(1):WNt_reward(2)),2)/twin1(2)-sum(r1(:,WNb_reward(1):WNb_reward(2)),2)/abs(bwin1(2)-bwin1(1));
    r2_ind = sum(r2(:,WNt_reward(1):WNt_reward(2)),2)/twin1(2)-sum(r2(:,WNb_reward(1):WNb_reward(2)),2)/abs(bwin1(2)-bwin1(1));
    
    cue_per_animal1(k) = mean(c1_ind); % mean cue and reward responses for each animal
    cue_per_animal2(k) = mean(c2_ind);
    reward_per_animal1(k) = mean(r1_ind);
    reward_per_animal2(k) = mean(r2_ind);
    
    filename_animal1 = ['cue_reward_correlation1_' animals{k}];
    filename_animal2 = ['cue_reward_correlation2_' animals{k}];
    
    % Trial-by-trial plot + linear regression for each animal
    regressionplot(c1_ind, r1_ind, 'Reward cue', 'Expected reward', filename_animal1, animals{k},  resdir)
    regressionplot(c2_ind, r2_ind, 'Punish cue', 'Surprising reward', filename_animal2, animals{k},  resdir)
end

regressionplot(cue_per_animal1', reward_per_animal1', 'Reward cue', 'Expected reward', 'cue_reward_correlation1', 'All animals', resdir) % all animals summarized
regressionplot(cue_per_animal2', reward_per_animal2', 'Punish cue', 'Surprising reward', 'cue_reward_correlation2', 'All animals', resdir)

% Regression plot for each cell individually
for c = 1:length(cellids)
    r1_cellid = BR_R1{c}; % get reward responses
    r2_cellid = BR_R2{c};
    
    c1_cellid = BR_C1{c}; % get cue responses
    c2_cellid = BR_C2{c};
    
    c1_full = sum(c1_cellid(:,WNt_cue(1):WNt_cue(2)),2) / abs(twin2(2)-twin2(1));  % non-corrected
    c1_bl = sum(c1_cellid(:,WNb_cue(1):WNb_cue(2)),2) / abs(bwin2(2)-bwin2(1));  % baseline
    c1_ind = (c1_full - c1_bl) ./ (c1_bl + 1);  % baseline-subtracted
    c2_full = sum(c2_cellid(:,WNt_cue(1):WNt_cue(2)),2) / abs(twin2(2)-twin2(1));  % non-corrected
    c2_bl = sum(c2_cellid(:,WNb_cue(1):WNb_cue(2)),2) / abs(bwin2(2)-bwin2(1));  % baseline
    c2_ind = (c2_full - c2_bl) ./ (c2_bl + 1);  % baseline-subtracted
    
    r1_full = sum(r1_cellid(:,WNt_reward(1):WNt_reward(2)),2) / abs(twin1(2)-twin1(1));  % non-corrected
    r1_bl = sum(r1_cellid(:,WNb_reward(1):WNb_reward(2)),2) / abs(bwin1(2)-bwin1(1));  % baseline
    r1_ind = (r1_full - r1_bl) ./ (r1_bl + 1);  % baseline-subtracted
    r2_full = sum(r2_cellid(:,WNt_reward(1):WNt_reward(2)),2) / abs(twin1(2)-twin1(1));  % non-corrected
    r2_bl = sum(r2_cellid(:,WNb_reward(1):WNb_reward(2)),2) / abs(bwin1(2)-bwin1(1));  % baseline
    r2_ind = (r2_full - r2_bl) ./ (r2_bl + 1);  % baseline-subtracted
    
    cellidt = cellids{c};
    cellidt(cellidt=='.') = '_';
    regressionplot(c1_ind, r1_ind, 'Reward cue', 'Expected reward', ['Cue_reward_correlation1_' cellidt], cellids{c},  resdir)
    regressionplot(c2_ind, r2_ind, 'Punish cue', 'Surprising reward', ['Cue_reward_correlation2_' cellidt], cellids{c},  resdir)
end

% -------------------------------------------------------------------------
function regressionplot(cue_resp, reward_resp, xtitle, ytitle, filenm, head, resdir)

% Plotting scatter plots
x = cue_resp;
y = reward_resp;

X = [ones(length(reward_resp),1) x];
[b,bint,r,rint,stats] = regress(y,X);
p = stats(3);
pR = corrcoef(y,X(:,2));
R = pR(3);
[b,stats] = robustfit(x,y);

% Regression plot
figure
scatter(x,y);
xlabel(xtitle)
ylabel(ytitle)
title(head)
% setmyplot_Balazs
axis square
icp = b(1);   % intercept
gr = b(2);   % gradient
xx = min(x):0.01:max(x);
yy = xx .* gr + icp;
hold on
plot(xx,yy,'Color',[0.6627 0.6196 0.4039],'LineWidth',2)   % overlay regression line
text('Units','normalized','Position',[0.7 0.7],...
    'String',{['p = ' num2str(p)] ['R = ' num2str(R)]})
set(gcf, 'renderer', 'painters')
savename = fullfile(resdir, [filenm '.fig']);
savename2 = fullfile(resdir, [filenm '.eps']);
savename3 = fullfile(resdir, [filenm '.jpg']);
saveas(gcf, savename)
saveas(gcf, savename2)
saveas(gcf, savename3)
close(gcf)