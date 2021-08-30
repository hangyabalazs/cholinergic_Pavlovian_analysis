function [stats] = response_correlation(cellid, var1, var2, resdir)
%   RESPONSE_CORRELATION(CELLID, VAR1, VAR2, RESDIR) response correlation
%   analysis.
%   RESPONSE_CORRELATION calculates correlation between responses defined by
%   VAR1 and VAR2 of a cell (CELLID). Result scatter plots are saved to
%   RESDIR. R^2 statistics, F statistic, p value and error variance are calculated as well (STATS).
%
%   See also ULTIMATE_PSTH and REGRESS

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   19-May-2020

narginchk(0,4)

if nargin < 4
    resdir = fullfile(getpref('cellbase', 'datapath'), '_paper_figs', 'response_correlation');
end

if ~isfolder(resdir) % make results directory
    mkdir(resdir);
end

% Time windows
wn = [-4 4];
dt = 0.001;  % time resolution
rspnswn_cue = [0 0.5];  % 'response window' - for calculating FR to be fitted with the TDRL model
rspnswn_cue_length = diff(rspnswn_cue);  % window length in s
rspnswn_reinforcement = [0 0.2];  % we allow different window after cue and reinforcement
rspnswn_reinforcement_length = diff(rspnswn_reinforcement);  % window length in s
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector in ms
testwin_cue_inx = find(time == rspnswn_cue(1)*1000):find(time == rspnswn_cue(2)*1000);   % the same for test windows
testwin_reinforcement_inx = find(time == rspnswn_reinforcement(1)*1000):find(time == rspnswn_reinforcement(2)*1000);

switch var1 %select input parameters for binraster calculation
    case 'cue'
        alevent1 = 'StimulusOn';
    case 'reward'
        alevent1 = 'DeliverAllFeedback';
        fevent1 = 'AllReward';
    case 'punish'
        alevent1 = 'DeliverAllFeedback';
        fevent1 = 'Punishment';
end

switch var2
    case 'cue'
        alevent2 = 'StimulusOn';
        fevent2 = fevent1;
    case 'reward'
        alevent2 = 'DeliverAllFeedback';
        fevent2 = 'AllReward';
        fevent1 = fevent2;
    case 'punish'
        alevent2 = 'DeliverAllFeedback';
        fevent2 = 'Punishment';
        fevent1 = fevent2;
end


if isequal(alevent1,'StimulusOn')  % use the approrpate response window
    testwin_inx1 = testwin_cue_inx;  % for cue
    rspnswn_length1 = rspnswn_cue_length;
elseif isequal(alevent1,'DeliverAllFeedback')
    testwin_inx1 = testwin_reinforcement_inx;  % for reinforcement (reward and punishment)
    rspnswn_length1 = rspnswn_reinforcement_length;
else
    error('MATLAB:hdb_pavlovian:cholinergic_RFmodel_main:aligneventNotRecognized',...
        'Aligning event not supported.')
end

if isequal(alevent2,'StimulusOn')  % use the approrpate response window
    testwin_inx2 = testwin_cue_inx;  % for cue
    rspnswn_length2 = rspnswn_cue_length;
elseif isequal(alevent2,'DeliverAllFeedback')
    testwin_inx2 = testwin_reinforcement_inx;  % for reinforcement (reward and punishment)
    rspnswn_length2 = rspnswn_reinforcement_length;
else
    error('MATLAB:hdb_pavlovian:cholinergic_RFmodel_main:aligneventNotRecognized',...
        'Aligning event not supported.')
end

% Calculate bin raster with 'ultimate_psth'
spike_data = cell(2,2);
for i = 1:2 %examine both trial types
    [~, ~, ~, ~, spt1, ~] = ultimate_psth(cellid, 'trial', alevent1, wn,...
        'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
        'event_filter','custom', 'filterinput', [fevent1 '==' num2str(i)],'maxtrialno',Inf,...
        'relative_threshold',0.1);
    
    [~, ~, ~, ~, spt2, ~] = ultimate_psth(cellid, 'trial', alevent2, wn,...
        'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',2,...
        'event_filter','custom', 'filterinput', [fevent2 '==' num2str(i)],'maxtrialno',Inf,...
        'relative_threshold',0.1);
    
    numTrials1 = size(spt1,1);  % number of trials
    [spike_count1] = nan(1,numTrials1);
    for j = 1:numTrials1
        spike_count1(j) = sum(spt1(j,(testwin_inx1))) / rspnswn_length1;   % firing rate in test window
    end
    spike_data{i,1} = spike_count1;   % firing rates
    
    
    numTrials2 = size(spt2,1);  % number of trials
    [spike_count2] = nan(1,numTrials2);
    for j = 1:numTrials2
        spike_count2(j) = sum(spt2(j,(testwin_inx2))) / rspnswn_length2;   % firing rate in test window
    end
    spike_data{i,2} = spike_count2;   % firing rates
    
    
    %Correlation for likely reward trials
    x=spike_data{i,1}';
    y = spike_data{i,2}';
    
    X = [ones(length(spike_data{i,1}),1) x];
    [b,bint,r,rint,stats] = regress(y,X);
    p = stats(3);
    pR = corrcoef(y,X(:,2));
    R = pR(3);
    [b,stats] = robustfit(x,y);
    % Regression plot
    figure
    scatter(x,y);
    xlabel([fevent1 'response'])
    ylabel('Cue response')
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
    fnm = [cellid '_' var1 '_' var2 '_correlation.fig'];
    fnm2 = [cellid '_' var1 '_' var2 '_correlation.jpg'];
    saveas(gcf, [resdir '\' fnm]);
    saveas(gcf, [resdir '\' fnm2]);
    close(gcf)
end

