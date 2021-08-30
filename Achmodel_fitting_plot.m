function model_values = Achmodel_fitting_plot(cellid,a1, c_a1, a2, c_a2, S, c_S, PARAM, c_PARAM, fr,err,control_err, baseline, resdir)
%ACHMODEL_FITTING_PLOT   Boxplot illustrating reinforcement model fitting on cholinergic neurons.
%   Y = ACHMODEL_FITTING_PLOT(CELLID,A1,A2,S,PARAM,FR,ERR,CONTROL_ERR,RESDIR)
%   creates boxplots illustrating the results of reinforcement learning
%   model fitting on cholinergic neuronal data. Four boxplots are plotted:
%   two for reward and punishment predicting cue and following
%   reinforcement with original contigencies of the behavioral task and
%   two with reversed contingencies. Median and average firing rate in the
%   response window, model results, FR in individual trials, median and
%   average firing rate of baseline firing are all indicated on the plots.
%       Input arguments:
%           CELLID - list of cell IDs to be plotted
%           A1, A2, S - best fit model parameters (for original task contingencies)
%           C_A1, C_A2, C_S - best fit model parameters (for reversed task contingencies)
%           PARAM, C_PARAM - task contingencies (original and reversed, respectively)
%           FR - firing rate in the response window
%           ERR, CONTROL_ERR - model fitting error (for original and
%           	reversed task contingencies)
%           BASELINE - baseline firing rate
%           RESDIR - results directory
%
%   Fitted values for cue, reward and punishment (Y, rows) are returned. 
%   Cells represent TrialType1 (likely reward) and TrialType2 (likely 
%   punishment) trials with original and reversed task contingencies.
%
%   See also ACH_RL_MODELFUN and BEST_FIT_CHOLINERGIC.

% Task contingencies
tE_R = {0.8, 0.25};  % task contingencies - reward probability in TrialType 1 and 2
tE_P = {0.1, 0.65};  % task contingencies - punishment probability in TrialType 1 and 2

% Parse trials
[inx, inx2] = getcontingencies(PARAM, {tE_R, tE_P}, false); % find different cues and outcomes in parameter matrix
[c_inx, c_inx2] = getcontingencies(c_PARAM, {tE_R, tE_P},true);

% Evaluate the model with the best fit parameters
modelvals1 = Ach_RL_modelfun([a1 a2 S], ...
    [tE_R{1} tE_P{1} 0 0 0; tE_R{1} tE_P{1} 1 0 0; tE_R{1} tE_P{1} 0 1 0]);  % fitted values for TrialType1 (likely reward)
modelvals2 = Ach_RL_modelfun([a1 a2 S], ...
    [tE_R{2} tE_P{2} 0 0 0; tE_R{2} tE_P{2} 1 0 0; tE_R{2} tE_P{2} 0 1 0]);  % fitted values for TrialType2 (likely punishment)

modelvals3 = Ach_RL_modelfun([c_a1 c_a2 c_S], ...
    [tE_P{1} tE_R{1} 0 0 0; tE_P{1} tE_R{1} 1 0 0; tE_P{1} tE_R{1} 0 1 0]);  % fitted values for TrialType1 (likely reward) - reversed contingencies
modelvals4 = Ach_RL_modelfun([c_a1 c_a2 c_S], ...
    [tE_P{2} tE_R{2} 0 0 0; tE_P{2} tE_R{2} 1 0 0; tE_P{2} tE_R{2} 0 1 0]);  % fitted values for TrialType2 (likely punishment) - reversed contingencies

model_values = {modelvals1 modelvals2 modelvals3 modelvals4}; % fitted values of the model (for T1 and T2 with original and reversed task contingencies)

% Plotting
figure
hold on
name = {'Cue' 'Reward' 'Punish'}; % x axis labels
subplot_title = {'Trial Type1' 'Trial Type2' 'Trial Type1 control' 'Trial Type2 control'};
color = {[0 0 1] [0 1 0] [1 0 0]}; % color code: blue- cue, green- reward, red- punishment
a = 0.1;
b = [1 2 3]; % vector helps plotting individual trials
ymax = max(fr); % min and max values for setting Y axis on subplots
ymin = min(fr);

for k = 1:4 % Plotting all subplots
    if mod(k,2)==1 % for subplot 1 and 3 (Trial Type 1)
        cinx = inx;
    else
        cinx = inx2; % for subplot 2 and 4 (Trial Type 2)
    end
    
    subplot(2,2,k) % Plot cue and reinforcements
%     boxplot([fr(cinx{1}); fr(cinx{2}); fr(cinx{3})],[zeros(size(fr(cinx{1}))); ones(size(fr(cinx{2}))); repmat(2,size(fr(cinx{3})))],'labels',name);
        bar([1 2 3], [mean(fr(cinx{1})) mean(fr(cinx{2})) mean(fr(cinx{3}))]);
    hold on
    errorbar([1 2 3], [mean(fr(cinx{1})) mean(fr(cinx{2})) mean(fr(cinx{3}))], [std(fr(cinx{1}))/sqrt(length(fr(cinx{1}))) std(fr(cinx{2}))/sqrt(length(fr(cinx{2}))) std(fr(cinx{3}))/sqrt(length(fr(cinx{3})))])
    hold on
    
    x1 = a .* randn(size(fr(cinx{1}),1),1) + b(1); % add a small random number to x value for each trial to avoid overlapping plotting
    x2 = a .* randn(size(fr(cinx{2}),1),1) + b(2);
    x3 = a .* randn(size(fr(cinx{3}),1),1) + b(3);
    xvalues = { x1 x2 x3 };
    for j = 1:3
        plot(xvalues{j}, fr(cinx{j}), 'o', 'MarkerSize',6,'MarkerEdgeColor', [color{j}]); % plot individual trials
        hold on
        plot(j, model_values{k}(j), 'o', 'MarkerSize',12,'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [color{j}]); % plot fitted model values
    end
    
    ylabel('Spike count');
    xlabel(subplot_title{k});
    set(gca,'YLim',[ymin-1 ymax+1]);
    set(gca, 'XLim', [0.8 3.2]);
end
% indicate fitting error (for both original and control contingencies)
% -placement should be improved
text(0, ymax+5, ['error: ' num2str(err) '     control err: ' num2str(control_err)], 'LineWidth', 2, 'Color','red');

% Save plot
maximize_figure;
cellid(cellid=='.') = '_';
saveas(gcf, [resdir '\' cellid '_spikenum_fitting_boxplot.fig']);
set(gcf, 'renderer', 'painters')
saveas(gcf, [resdir '\' cellid '_spikenum_fitting_boxplot.eps']);
saveas(gcf, [resdir '\' cellid '_spikenum_fitting_boxplot.jpeg']);
close(gcf);

% -------------------------------------------------------------------------
function [inx, inx2] = getcontingencies(PARAM, cont, isreversed)

if ~isreversed % if task contingencies are reversed (control)
    E_R = cont{1};
else
    E_R = cont{2};
end

rewardcue = PARAM(:,3)==0 & PARAM(:,4)==0 & PARAM(:,1)==E_R{1};   % likely reward cues
punishcue = PARAM(:,3)==0 & PARAM(:,4)==0 & PARAM(:,1)==E_R{2};   % likely punishment cues
ereward= PARAM(:,3)==1 & PARAM(:,4)==0 & PARAM(:,1)==E_R{1};   % expected reward
spunish= PARAM(:,3)==0 & PARAM(:,4)==1 & PARAM(:,1)==E_R{1};   % surprising punishment
sreward= PARAM(:,3)==1 & PARAM(:,4)==0 & PARAM(:,1)==E_R{2};   % surprising reward
epunish= PARAM(:,3)==0 & PARAM(:,4)==1 & PARAM(:,1)==E_R{2};   % surprising punishment

inx = {rewardcue ereward spunish};
inx2 = {punishcue sreward epunish};