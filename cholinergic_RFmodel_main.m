function [a1, a2, S, err, c_a1, c_a2, c_S, c_err, model_values] = ...
    cholinergic_RFmodel_main(cellids, window, varargin)
%CHOLINERGIC_RFMODEL_MAIN   Fitting a reinforcement learning model to
%electrophysiological recordings of cholinergic basal forebrain neurons.
%   [A1, A2, S, ERR] = CHOLINERGIC_RFMODEL_MAIN(CELLIDS, WINDOW, VARARGIN)
%   is the main function for reinforcement learning model fitting on
%   cholinergic neuronal activity. It creates a binraster for each neuron
%   and fits the RF model on the data. A1, A2 and S best fit paramaters of
%   the model are returned as output. ERR is the MLE error of the best fit
%   model. Eventually, it creates an average boxplot and CDF for neural responses
%   and fitted/simulated values
%
%   Input arguments:
%       CELLIDS - list of cell IDs to be analyzed
%       WINDOW - full data windoow around the aligning event included in
%           the analysis
%
%   Optional input arguments:
%       DT - deafult, 0.001; temporal resolution for bin raster calculation
%           in seconds
%       BASELINEWIN - default, [-3 -2]; time window for baseline firing
%           rate relative to the aligning event, in seconds; baseline
%           window is only used for plotting and not for the model fitting
%       RESPONSEWIN_CUE - default, [0 0.5]; time window for firing rate
%           after cue presentation, in seconds
%       RESPONSEWIN_REINFOCEMENT - default, [0 0.2]; time window for firing
%           rate after reinforcemnt presentation (reward or punishment), in
%           seconds
%       RESDIR - default, current working directory; results directory for
%           saving
%
%   See also ULTIMATE_PSTH and BEST_FIT_CHOLINERGIC

%   Panna Hegedus, Balazs Hangya
%   Insititute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   05-Feb-2020

%   Code review: BH 4/16/20

% Parsing inputs
prs = inputParser;
addRequired(prs,'cellids',@(s)isempty(s)|iscell(s)|iscellstr(s)|ischar(s)) % cellIDs - cells
addRequired(prs,'window',@(s)isnumeric(s)&isequal(length(s),2))  % full data window in seconds
addParameter(prs,'dt',0.001,@isnumeric)   % resolution, in seconds
addParameter(prs,'baselinewin',[-3 -2],@(s)isnumeric(s)&isequal(length(s),2))   % baseline window in seconds
addParameter(prs,'responsewin_cue',[0 0.5],@(s)isnumeric(s)&isequal(length(s),2))   % response window for cue, in seconds
addParameter(prs,'responsewin_reinforcement',[0 0.2],@(s)isnumeric(s)&isequal(length(s),2))   % response window for reinforcement, in seconds
addParameter(prs,'resdir',pwd,@(s)ischar(s))   % results directory
parse(prs,cellids,window,varargin{:})
g = prs.Results;

% Cells
if isempty(g.cellids) % if no cellIDs were provided
    %     choosecb('Cholinergic_pavlovian')
    selstr = '"ChAT+"==1';
    cellids = selectcell(selstr);   % select cholinergic neurons
end

if ~isfolder(g.resdir)
    mkdir(g.resdir)
end

% Time windows
wn = g.window;  % full data window
dt = g.dt;  % time resolution
bslnwn = g.baselinewin;   % basseline window for model fitting (only for plotting purposes)
bsl_length = diff(g.baselinewin);   % baseline window length in s
rspnswn_cue = g.responsewin_cue;  % 'response window' - for calculating FR to be fitted with the TDRL model
rspnswn_cue_length = diff(g.responsewin_cue);  % window length in s
rspnswn_reinforcement = g.responsewin_reinforcement;  % we allow different window after cue and reinforcement
rspnswn_reinforcement_length = diff(g.responsewin_reinforcement);  % window length in s
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector in ms
baseline_inx = find(time == bslnwn(1)*1000):find(time == bslnwn(2)*1000);   % baseline window indices wrt. to full data window
testwin_cue_inx = find(time == rspnswn_cue(1)*1000):find(time == rspnswn_cue(2)*1000);   % the same for test windows
testwin_reinforcement_inx = find(time == rspnswn_reinforcement(1)*1000):find(time == rspnswn_reinforcement(2)*1000);

% Calculate firing rates
NumCells = length(cellids);   % number of analyzed cells
[a1, a2, S, err] = deal(nan(1,NumCells));
[c_a1, c_a2, c_S, c_err] = deal(nan(1,NumCells));
[BasFR, RespFR, maxval_all] = deal(cell(1,NumCells));
[T1_spike_avg, T2_spike_avg] = deal(nan(NumCells, 3));
model_values = cell(1,NumCells);

for i = 1:NumCells   % loop through cells
    cellid = cellids{i}; % current cell
    
    eps = zeros(1,3);
    [spike_data, baseline_data, maxvalue] = deal(cell(1,3));
    alignevent = {'StimulusOn', 'DeliverAllFeedback', 'DeliverAllFeedback'};
    filterevent = {'StimulusOn', 'AllReward==1 | AllReward==2', 'Punishment==1 | Punishment==2'};
    
    for k = 1:3  % loop through cue, reward and punishment
        alevent = alignevent{k};
        fevent = filterevent{k};
        if isequal(alevent,'StimulusOn')  % use the approrpate response window
            testwin_inx = testwin_cue_inx;  % for cue
            rspnswn_length = rspnswn_cue_length;
        elseif isequal(alevent,'DeliverAllFeedback')
            testwin_inx = testwin_reinforcement_inx;  % for reinforcement (reward and punishment)
            rspnswn_length = rspnswn_reinforcement_length;
        else
            error('MATLAB:hdb_pavlovian:cholinergic_RFmodel_main:aligneventNotRecognized',...
                'Aligning event not supported.')
        end
        
        % Calculate bin raster with 'ultimate_psth'
        [~, ~, ~, ~, spt, stats] = ultimate_psth( cellid, 'trial', alevent, wn,...
            'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
            'event_filter','custom', 'filterinput', fevent,'maxtrialno',Inf,...
            'relative_threshold',0.1);
        
        numTrials = size(spt,1);  % number of trials
        [spike_count, baseline] = deal(nan(1,numTrials));
        for j = 1:numTrials
            baseline(j) = sum(spt(j, baseline_inx)) / bsl_length;   % firing rate in baseline window
            spike_count(j) = sum(spt(j,(testwin_inx))) / rspnswn_length;   % firing rate in test window
        end
        eps(k) = nanstd(spike_count);   % SD of firing rate (for the noise model in MLE fitting)
        spike_data{k} = spike_count;   % firing rates
        baseline_data{k} = baseline;   % baseline firing rates
        maxvalue(k) = stats.maxvalue;
    end
    maxval_all{i} = maxvalue;
    % Fit TDRL model
    [a1(i), a2(i), S(i), PARAM, fr, err(i)] = best_fit_cholinergic(cellid,spike_data,eps);
    
    % Fit TDRL model with reverse trial continencies
    [c_a1(i), c_a2(i), c_S(i), c_PARAM, c_fr, c_err(i)] = ...
        best_fit_cholinergic(cellid,spike_data,eps,'isreversed',true);
    
    % Plot model fitting results
    sim_values = Achmodel_fitting_plot(cellid,a1(i), c_a1(i), a2(i), c_a2(i), S(i), c_S(i), PARAM, c_PARAM,fr,err(i),c_err(i), baseline ,g.resdir);
    model_values{i} = horzcat(sim_values{:}); % reorganize simulated values
    
    % Trial events
    VE = loadcb(cellid,'TrialEvents');   % load events
    
    % Make indices for reward and punishment trials
    VE.AllReward(any(isnan(VE.AllReward),1))=[];
    reward_inx = VE.AllReward;
    VE.Punishment(any(isnan(VE.Punishment),1))=[];
    punish_inx = VE.Punishment;
    
    % Average spike counts
    T1_spike_avg(i,1)= mean(spike_data{1}(VE.TrialType==1));
    T1_spike_avg(i,2)= mean(spike_data{2}(reward_inx==1));
    T1_spike_avg(i,3)= mean(spike_data{3}(punish_inx==1));
    
    T2_spike_avg(i,1)= mean(spike_data{1}(VE.TrialType==2));
    T2_spike_avg(i,2)= mean(spike_data{2}(reward_inx==2));
    T2_spike_avg(i,3)= mean(spike_data{3}(punish_inx==2));
end

%%%%% PLOT AVERAGE CHOLINERGIC CELL RESPONSE %%%%%%%%%

% Average model parameters and simulated values
a1_avg = mean(a1);
a2_avg = mean(a2);
S_avg = mean(S);

c_a1_avg = mean(c_a1);
c_a2_avg = mean(c_a2);
c_S_avg = mean(c_S);

% Simulate average firing rate
% Task contingencies
tE_R = {0.8, 0.25};  % task contingencies - reward probability in TrialType 1 and 2
tE_P = {0.1, 0.65};  % task contingencies - punishment probability in TrialType 1 and 2

% Evaluate the model with the average best fit parameters
modelvals1 = Ach_RL_modelfun([a1_avg a2_avg S_avg], ...
    [tE_R{1} tE_P{1} 0 0 0; tE_R{1} tE_P{1} 1 0 0; tE_R{1} tE_P{1} 0 1 0]);  % fitted values for TrialType1 (likely reward)
modelvals2 = Ach_RL_modelfun([a1_avg a2_avg S_avg], ...
    [tE_R{2} tE_P{2} 0 0 0; tE_R{2} tE_P{2} 1 0 0; tE_R{2} tE_P{2} 0 1 0]);  % fitted values for TrialType2 (likely punishment)

modelvals3 = Ach_RL_modelfun([c_a1_avg c_a2_avg c_S_avg], ...
    [tE_P{1} tE_R{1} 0 0 0; tE_P{1} tE_R{1} 1 0 0; tE_P{1} tE_R{1} 0 1 0]);  % fitted values for TrialType1 (likely reward) - reversed contingencies
modelvals4 = Ach_RL_modelfun([c_a1_avg c_a2_avg c_S_avg], ...
    [tE_P{2} tE_R{2} 0 0 0; tE_P{2} tE_R{2} 1 0 0; tE_P{2} tE_R{2} 0 1 0]);  % fitted values for TrialType2 (likely punishment) - reversed contingencies

avg_model_values = {modelvals1 modelvals2 modelvals3 modelvals4}; % fitted values of the model (for T1 and T2 with original and reversed task contingencies)

% Average plot
x1 = 0.1 .* randn(size(T1_spike_avg,1),1) + 1; % add a small random number to x value for each trial to avoid overlapping plotting
x2 = 0.1 .* randn(size(T1_spike_avg,1),1) + 2;
x3 = 0.1 .* randn(size(T1_spike_avg,1),1) + 3;
xvalues = { x1 x2 x3 };

color = {[0 0 1] [0 1 0] [1 0 0]}; % color code: blue- cue, green- reward, red- punishment

for a = 1:4
    if mod(a,2)==0
        spike_avg = T2_spike_avg;
    else
        spike_avg = T1_spike_avg;
    end
    
    subplot(2,2,a)
%     boxplot([spike_avg(:,1); spike_avg(:,2); spike_avg(:,3)],[zeros(size(spike_avg(:,1))); ones(size(spike_avg(:,2))); repmat(2,size(spike_avg(:,3)))],'labels',{'Cue' 'Reward' 'Punishment'})
      bar([1 2 3], [mean(spike_avg(:,1)) mean(spike_avg(:,2)) mean(spike_avg(:,3))])
    hold on
          errorbar([1 2 3], [mean(spike_avg(:,1)) mean(spike_avg(:,2)) mean(spike_avg(:,3))],[std(spike_avg(:,1))/sqrt(length(spike_avg(:,1))) std(spike_avg(:,2))/sqrt(length(spike_avg(:,2))) std(spike_avg(:,3))/sqrt(length(spike_avg(:,3)))]);
hold on
    for x = 1:3
        plot(xvalues{x}, spike_avg(:,x), 'o', 'MarkerSize',6,'MarkerEdgeColor', [color{x}]); % plot individual cells
        plot(x, mean(spike_avg(:,x)), '*', 'MarkerSize',12,'LineWidth', 2, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]); % plot average neuronal response
        plot(x, avg_model_values{a}(x), 'o', 'MarkerSize',12,'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [color{x}]); % plot fitted model values
    end
end
set(gcf, 'renderer', 'painters')
saveas(gcf, fullfile(g.resdir, 'average_responses.eps')) % save plot
saveas(gcf, fullfile(g.resdir, 'average_responses.jpg')) % save
close(gcf)

% Check model parameter consistency within animals
animals = getvalue('RatID_tag',cellids);
animalpairs = [];   % indices for pairs of cells recorded in the same mouse
nonanimalpairs = [];   % indices for pairs of cells recorded in different mice
for iA = 1:NumCells
    mtch = find(strcmp(animals{iA},animals(iA+1:end)));
    if ~isempty(mtch)
        animalpairs = [animalpairs; repmat(iA,length(mtch),1) iA+mtch]; %#ok<AGROW>
    end
    nonanimalpairs = [nonanimalpairs; repmat(iA,NumCells-iA-length(mtch),1) setdiff((iA+1:NumCells)',iA+mtch)]; %#ok<AGROW>
end
a1ap = a1(animalpairs);   % eta1 parameters for animal-pair cells
a1anp = a1(nonanimalpairs);   % eta1 parameters for non-animal-pair cells
a2ap = a2(animalpairs);   % eta2 parameters for animal-pair cells
a2anp = a2(nonanimalpairs);   % eta2 parameters for nonanimal-pair cells
boxstat(abs(diff(a1ap,[],2)),abs(diff(a1anp,[],2)),'animal-pairs','non-animal-pairs')  % compare pairwise differences in eta1 within and across animals
title('Compare pairwise differences in eta1 within and across animals')
set(gcf, 'renderer', 'painters')
saveas(gcf, fullfile(g.resdir, 'within_animal_differences_eta1.eps')) % save plot
saveas(gcf, fullfile(g.resdir, 'within_animal_differences_eta1.jpg')) % save
saveas(gcf, fullfile(g.resdir, 'within_animal_differences_eta1.fig')) % save
close(gcf)
boxstat(abs(diff(a2ap,[],2)),abs(diff(a2anp,[],2)),'animal-pairs','non-animal-pairs')  % compare pairwise differences in eta2 within and across animals
title('Compare pairwise differences in eta2 within and across animals')
set(gcf, 'renderer', 'painters')
saveas(gcf, fullfile(g.resdir, 'within_animal_differences_eta2.eps')) % save plot
saveas(gcf, fullfile(g.resdir, 'within_animal_differences_eta2.jpg')) % save
saveas(gcf, fullfile(g.resdir, 'within_animal_differences_eta2.fig')) % save
close(gcf)

% Check model parameter consistency within sessions
sessions = getvalue('SessionID_tag',cellids);
sessionpairs = [];   % indices for pairs of cells recorded in the same session
nonsessionpairs = [];   % indices for pairs of cells recorded in different sessions
for iS = 1:NumCells
    mtch = find(strcmp(sessions{iS},sessions(iS+1:end)));
    if ~isempty(mtch)
        sessionpairs = [sessionpairs; repmat(iS,length(mtch),1) iS+mtch]; %#ok<AGROW>
    end
    nonsessionpairs = [nonsessionpairs; repmat(iS,NumCells-iS-length(mtch),1) setdiff((iS+1:NumCells)',iS+mtch)]; %#ok<AGROW>
end
a1p = a1(sessionpairs);   % eta1 parameters for session-pair cells
a1np = a1(nonsessionpairs);   % eta1 parameters for non-session-pair cells
a2p = a2(sessionpairs);   % eta2 parameters for session-pair cells
a2np = a2(nonsessionpairs);   % eta2 parameters for nonsession-pair cells
boxstat(abs(diff(a1p,[],2)),abs(diff(a1np,[],2)),'session-pairs','non-session-pairs')  % compare pairwise differences in eta1 within and across sessions
title('Compare pairwise differences in eta1 within and across sessions')
set(gcf, 'renderer', 'painters')
saveas(gcf, fullfile(g.resdir, 'within_session_differences_eta1.eps')) % save plot
saveas(gcf, fullfile(g.resdir, 'within_session_differences_eta1.jpg')) % save
saveas(gcf, fullfile(g.resdir, 'within_session_differences_eta1.fig')) % save
close(gcf)
boxstat(abs(diff(a2p,[],2)),abs(diff(a2np,[],2)),'session-pairs','non-session-pairs')  % compare pairwise differences in eta2 within and across sessions
title('Compare pairwise differences in eta2 within and across sessions')
set(gcf, 'renderer', 'painters')
saveas(gcf, fullfile(g.resdir, 'within_session_differences_eta2.eps')) % save plot
saveas(gcf, fullfile(g.resdir, 'within_session_differences_eta2.jpg')) % save
saveas(gcf, fullfile(g.resdir, 'within_session_differences_eta2.fig')) % save
close(gcf)

% CDF of simulated and recorded neural responses
cholinergic_CDF(T1_spike_avg, T2_spike_avg, model_values, g.resdir)

% Save workspace variables
save(fullfile(g.resdir,'TDRLmodel_output_variables.mat'),'a1','a2','S','err',...
    'c_a1','c_a2','c_S','c_err','c_fr','BasFR','RespFR','model_values');