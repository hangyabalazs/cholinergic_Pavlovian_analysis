function cholinergic_analysis_pavlovian_main(achcells, preproc, fig1spec, fig2spec, fig3spec, fig4spec, extraspec)
%CHOLINERGIC_ANALYSIS_PAVLOVIAN_MAIN(ACHCELLS, PREPROC, FIG1SPEC) Main analysis function for
%cholinergic pavlovian project.
%   CHOLINERGIC_ANALYSIS_PAVLOVIAN_MAIN(ACHCELLS, PREPROC, FIG1SPEC) is the
%   main function for basal forebarin cholinergic neural data analysis.
%   Data analysis is performed for CellIDs in ACHCELLS - or all
%   well-isolated units determined by SELECT_ACH_CELLS if the first input
%   is empty. The analysis modules are controlled by logical input
%   variables grouped by figures (and pre-processing including
%   auto-correlation and cross-correlation and analysis aiding the decision
%   about optogenetic tagging). This code is for further analysis after
%   running QUICKANALYSIS_PAVLOVIAN2_P on the CellBase dataset.
%   Input parameters:
%       PREPROC - 1x3 logical array controlling data preprocessing
%           (autocorrelation, cross-correlation and statistics on optical
%           tagging efficiency, respectively).
%       FIG1SPEC - 1x3 logical array controlling analysis for behavior,
%           optogenetic tagging and cluster quality.
%       FIG2SPEC - 1x2 logical array controlling analysis for neuronal response to
%           cue, reward and punishment.
%       FIG3SPEC - 1x1 logical input controlling the fitting of a reinforcement model
%           to neuronal data.
%       FIG4SPEC - 1x4 logical input controlling the analysis for cue,
%           anticipatory lick response and reaction time
%       EXTRASPEC - 1x4 logical input controlling extra analyses inclusing
%           CCG pairs, PSTH by area and PSTH conditioned on outcome history
%           and cue response size
%
%   See also QUICKANALYSIS_PAVLOVIAN2_P, SELECT_ACH_CELLS, TAGGEDPROP, ACH_ACG
%   and ACH_CCG.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   29-Apr-2020

% Input argument check
narginchk(0,6)
if nargin < 1
    achcells = [];
end

% Choose CellBase
usr = getenv('username');
if ismember(usr,{'hangya.balazs','hangyab','Hangya Balázs'})
    choosecb('Pavlovian_cholinergic')   % choose CellBase
else
    choosecb('Cholinergic_pavlovian')   % choose CellBase
end

% Control of different analysis domains
if nargin < 2
    preproc = [1 1 1];
end
if nargin < 3
    fig1spec = [1 1 1];
end
if nargin < 4
    fig2spec = [1 1];
end
if nargin < 5
    fig3spec = 1;
end
if nargin < 6
    fig4spec = [1 1 1 1];
end
if nargin < 7
    extraspec = [0 0 0 0];
end

% Parse inputs
perform_acg = preproc(1);   % ACG
perform_ccg = preproc(2);   % CCG
tagging_stat = preproc(3);   % TAGGING STATISTICS

behav_analysis = fig1spec(1);
opto_tagging = fig1spec(2);
clustering = fig1spec(3);

example_PSTH_plot = fig2spec(1);
avg_PSTH_plot = fig2spec(2);

cue_rt = fig4spec(1);
cue_lick_corr = fig4spec(2);
islickPSTH = fig4spec(3);
rt_lick_corr = fig4spec(4);

psth_by_area = extraspec(1);
ccg_pairs = extraspec(2);
cue_reward_corr = extraspec(3);
more_PSTHs = extraspec(4);

% Getting cellids of cholinergic BF neurons
loadcb
achcells = select_ach_cells(achcells);
location = getvalue('Area1', achcells);
achcells = setdiff(achcells, achcells(strcmp(location, 'AcbC'))) % remove the cholinerg cell from accumbens

% Results directory
resdir = [getpref('cellbase','datapath') '\_paper_figs_2021MARCH\']; % directory for saving figures/analysis results
if ~isfolder(resdir)
    mkdir(resdir)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate autocorrelograms
acg_resdir = fullfile(resdir,'acg');   % results directory
if perform_acg
    ach_acg(achcells, acg_resdir,true);
end

% Calculate cross-correlograms
ccg_resdir = fullfile(resdir,'ccg');   % results directory
if perform_ccg
    ach_ccg(achcells, ccg_resdir,true);
end

% Perform statistics to decide the efficacy of optogenetic tagging
if tagging_stat
    error_list = taggedprop(achcells, true);
    close all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT READY YET
% FIGURE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if behav_analysis
    
    % Example lick raster and PETH - FIG1.D
    resdir_lick_example = fullfile(resdir,'Fig1','lick_example');   % results directoryfigure;  % example lick PSTH + raster
    if ~isfolder(resdir_lick_example)
        mkdir(resdir_lick_example)
    end
    viewlick({'HDB36' '190518a'},'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav', 'sigma', 0.05,...
        'ShowEvents',{{'StimulusOn' 'StimulusOff'}},...  %DOUBLE CHECK SESSION
        'Partitions','#TrialType','window',[-3 3]) % lick rasters ans PSTH aligned to stimulus onset
    H = gcf;
    fnm = fullfile(resdir_lick_example, 'lick_example.fig');
    set(H, 'renderer', 'painters')
    fnm2 = fullfile(resdir_lick_example, 'lick_example.eps');
    fnm3 = fullfile(resdir_lick_example, 'lick_example.jpg');
    saveas(H,fnm);
    saveas(H,fnm2);
    saveas(H,fnm3);
    
    close(H)
    
    % Average lick PETH - FIG1.E-F
    resdir_lick_psth = fullfile(resdir,'Fig1','lick_psth');   % results directory
    Ach_lick_psth_summary_allcond(achcells,resdir_lick_psth,true)  % average lick PSTH
    Ach_lick_psth_summary(achcells,resdir_lick_psth,true)  % average lick PSTH
    Ach_lick_poppsth(achcells,resdir_lick_psth,true) % color-coded lick raster
    
    % Behavioral statistics - Wilcoxon signed rank test
    %     resdir_astat = fullfile(resdir,'Fig1','anticipatory_stat');
    %     T18 = anticipatory_stat('HDB18',resdir_astat);
    %     T25 = anticipatory_stat('HDB25',resdir_astat);
    %     T36 = anticipatory_stat('HDB36',resdir_astat);
    %     T38 = anticipatory_stat('HDB38',resdir_astat);
    %     T01 = anticipatory_stat('HW001',resdir_astat);
    %     behavioral_data = [{T18} {T25} {T36} {T38} {T01}];
    %     labels = {'HDB18' 'HDB25' 'HDB36' 'HDB38' 'HW001'};
    %     lineplot_with_errorbars(behavioral_data, labels, resdir_astat); % plot results: line plot with error bars for each animal
    
    % average spike shape
    resdir_spikeshape = fullfile(resdir,'Fig1','spike_shape_analysis');   % results directory
    issave = true;
    SpikeShape=spikeshapeanalysis_p(achcells,issave, resdir_spikeshape)
    
    avg_waveform = nan(length(achcells), 27);
    for a = 1:length(achcells)
        avg_waveform(a,:)=SpikeShape{a}.Spike
    end
    %plot avg waveform
    H = figure;
    errorshade([1:size(avg_waveform,2)],mean(avg_waveform),std(avg_waveform),...
        'LineColor',[0 0 1],'ShadeColor',[0 0 1])
    set(gcf, 'renderer', 'painters'); % set renderer
    
    % save
    saveas(gcf, fullfile(resdir_spikeshape, ['avg_spikeshape.eps']));
    saveas(gcf, fullfile(resdir_spikeshape, ['avg_spikeshape.jpg']));
end

if opto_tagging
    if ~tagging_stat % if statistics on light activation was not performed previously
        error_list = taggedprop(achcells, true);  % RESULTS DIRECTORY
    end
    resdir_LS = fullfile(resdir,'Fig1','light_stimulation'); % define result's directories
    LS_stat_dir = fullfile(getpref('cellbase','datapath'),'taggingsummary');
    
    % Example light evoked PSTH - FIG.1K
    cholinergic_LS({'HDB25_180526a_5.1'}, resdir_LS)
    
    % Population PSTH for light evoked PSTH - FIG.1K
    ChAT_lightpoppsth(achcells, LS_stat_dir, resdir_LS)
    
    % Latency and jitter CDF for light activation - FIG.1L
    ChAT_R_L_J(achcells, resdir_LS, 'tag', true)
end

if clustering
    
    resdir_C = fullfile(resdir,'Fig1','clustering'); % define result's directories
    
    % Cluster quality measures
    cluster_quality_stats(achcells, resdir_C, 'cholinergic')  % cholinergic
    
    sc = selectcell('"ID_PC">20&"Lr_PC"<0.15');  % cells with criterion quality
    ntcells = setdiff(sc,achcells);  % untagged population
    cluster_quality_stats(ntcells, resdir_C, 'untagged')  % untagged
    
    % Spike shape correlation
    NumCells = length(achcells);
    R = nan(1,NumCells);
    for iC = 1:NumCells
        R(iC) = spikeshapecorr(achcells{iC});  % should be fixed - nlx-based
    end
    median_R = median(R);
    SE_R = se_of_median(R);
    range_R = [min(R) max(R)];
    save(fullfile(resdir_C,'spikeshapecorr_variables.mat'),'median_R','SE_R',...
        'range_R','R');
    
    % H index distribution
    H = Hindex_distribution(ntcells,achcells);
    fnm = fullfile(resdir_C, 'tagging_histogram.fig');
    saveas(H,fnm);
    close(H)
    
    % Example cholinergic cluster
    [selcluster, neighborclusters, selLS] = example_cluster('HDB36_190512a_3.2');  % best example: HDB36_190512a_3.2; second best: HDB_36_190515a_3.1
    fnm = fullfile(resdir_C, 'example_cluster.mat');
    save(fnm,'selcluster','neighborclusters','selLS');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preselect neurons for further analysis
areas = getvalue('Area1', achcells); % select HDB neurons
hdbcells = achcells(strcmp(areas, 'HDB'));

% Filter: exclude cells with less then 10 surprising reward trials
NumCells = length(hdbcells);
delstr = ones(1,NumCells);
for i = 1:NumCells
    cellid = hdbcells(i);
    E = loadcb(cellid,'TrialEvents');   % load events
    if (sum(E.Reward==2) < 10) % less then 10 surprising reward trials
        delstr(i) = 0;
    end
end
selected_hdbcells = hdbcells(delstr==1);

% Filter: exclude cells with less then 10 surprising reward trials
NumCells = length(achcells);
delstr = ones(1,NumCells);
for i = 1:NumCells
    cellid = achcells(i);
    E = loadcb(cellid,'TrialEvents');   % load events
    if (sum(E.Reward==2) < 10) % less then 10 surprising reward trials
        delstr(i) = 0;
    end
end
selected_achcells = achcells(delstr==1);

delstr2 = ones(1,NumCells);
for i = 1:NumCells
    cellid = achcells(i);
    E = loadcb(cellid,'TrialEvents');   % load events
    if (sum(E.Omission==2) < 5) % less then 10 surprising reward trials
        delstr2(i) = 0;
    end
end
omission_achcells = achcells(delstr2==1);

% Example PSTHs
if example_PSTH_plot
    resdir_example=fullfile(resdir, 'Fig2', 'example');
    if ~isfolder(resdir_example)
        mkdir(resdir_example);
    end
    %     example_cells = {'HDB36_190426a_3.1', 'HDB36_190504a_5.1', 'HDB36_190508a_3.2'};
    example_cells = achcells;
    for e=1:length(example_cells)
        e_cell = example_cells{e};
        cellidt = e_cell;
        cellidt(cellidt=='.') = '_';
        G = figure
        viewcell2b(e_cell,'TriggerName','StimulusOn','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType','window',[-3 3])
        saveas(G, fullfile(resdir_example, ['cue_example_' cellidt '.jpg']))
        set(G, 'renderer', 'painters')
        saveas(G, fullfile(resdir_example, ['cue_example_' cellidt '.eps']))
        close(G)
        %
        %         H = figure
        %         viewcell2b(e_cell,'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#AllReward','window',[-3 3])
        %         saveas(H, fullfile(resdir_example, ['reward_example_' cellidt '.jpg']))
        %         set(H, 'renderer', 'painters')
        %         saveas(H, fullfile(resdir_example, ['reward_example_' cellidt '.eps']))
        %         close(H)
        %
        %
        %         I = figure
        %         viewcell2b(e_cell,'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#Punishment','window',[-3 3])
        %         saveas(I, fullfile(resdir_example, ['punish_example_' cellidt '.jpg']))
        %         set(I, 'renderer', 'painters')
        %         saveas(I, fullfile(resdir_example, ['punish_example_' cellidt '.eps']))
        %         close(I)
        J = figure
        viewcell2b(e_cell,'TriggerName','StimulusOn','SortEvent','TrialStart','sigma', 0.2,'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#Omission','window',[-3 3])
        saveas(J, fullfile(resdir_example, ['omission_example_cuealigned_' cellidt '.jpg']))
        set(J, 'renderer', 'painters')
        saveas(J, fullfile(resdir_example, ['omission_example_cuealigned_' cellidt '.eps']))
        close(J)
        
        %
        J = figure
        viewcell2b(e_cell,'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.2,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#Omission','window',[-3 3])
        saveas(J, fullfile(resdir_example, ['omission_example_' cellidt '.jpg']))
        set(J, 'renderer', 'painters')
        saveas(J, fullfile(resdir_example, ['omission_example_' cellidt '.eps']))
        close(J)
    end
    %later - HDB36_190426_3_1, HDB36_190504a_5_1, HDB36_190508a_3_2
end

% Average PSTHs
if avg_PSTH_plot % average PSTH, boxplot for response magnitude and maxvalue
    
    % Include all trials
    % All cholinergic neurons
    PSTHdir1 = 'avg_PSTH_nonadaptive2';
    PSTHdir1_selected = 'avg_PSTH_nonadaptive2_selected';
    [psth_R1, psth_R2, time] = avg_psth_cholinergic(achcells, 'cue', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'all_cholinergic'));
    avg_psth_cholinergic(achcells, 'reward', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'all_cholinergic'));
    avg_psth_cholinergic(achcells, 'punish', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'all_cholinergic'));
    avg_psth_cholinergic(achcells, 'omission', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'all_cholinergic'));
    avg_psth_cholinergic(selected_achcells, 'omission', 'none', fullfile(resdir, 'Fig2', PSTHdir1_selected, 'all_cholinergic'));
    avg_psth_cholinergic(omission_achcells, 'omission', 'none', fullfile(resdir, 'Fig2', PSTHdir1_selected, 'all_cholinergic'));
    
    
    label = getvalue('RatID_Tag',achcells);  % FIG.SX
    resdir_g1 = fullfile(resdir,'Fig2','group_by_animal','TrialType1');
    psth_by_label(psth_R1,time,label,resdir_g1);
    resdir_g1 = fullfile(resdir,'Fig2','group_by_animal','TrialType2');
    psth_by_label(psth_R2,time,label,resdir_g1);
    label = getvalue('Area1',achcells);
    resdir_g2 = fullfile(resdir,'Fig2','group_by_area','TrialType1');
    psth_by_label(psth_R1,time,label,resdir_g2);
    resdir_g2 = fullfile(resdir,'Fig2','group_by_area','TrialType2');
    psth_by_label(psth_R2,time,label,resdir_g2);
    
    % All HDB neurons
    avg_psth_cholinergic(hdbcells, 'cue', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'all_hdb'));
    avg_psth_cholinergic(hdbcells, 'reward', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'all_hdb'));
    avg_psth_cholinergic(hdbcells, 'punish', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'all_hdb'));
    
    % Filtered cholinergic neurons
    avg_psth_cholinergic(selected_achcells, 'cue', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_cholinergic'));
    [psth_R1, psth_R2, time] = avg_psth_cholinergic(selected_achcells, 'reward', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_cholinergic'));
    avg_psth_cholinergic(selected_achcells, 'punish', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_cholinergic'));
    
    %Quantify effect size
    effect_size(selected_achcells, 'cue', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_cholinergic_effect_size'));
    effect_size(selected_achcells, 'reward', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_cholinergic_effect_size'));
    effect_size(selected_achcells, 'punish', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_cholinergic_effect_size'));
    
    ChAT_R_L_J(achcells, fullfile(resdir, 'Fig2', PSTHdir1, 'all_cholinergic'), 'cue', true);
    ChAT_R_L_J(achcells, fullfile(resdir, 'Fig2', PSTHdir1, 'all_cholinergic'), 'rew', true);
    ChAT_R_L_J(achcells, fullfile(resdir, 'Fig2', PSTHdir1, 'all_cholinergic'), 'pun', true);
    
    ChAT_R_L_J(selected_achcells, fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_cholinergic'), 'cue', true);
    ChAT_R_L_J(selected_achcells, fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_cholinergic'), 'rew', true);
    ChAT_R_L_J(selected_achcells, fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_cholinergic'), 'pun', true);
    psth_heatmap(psth_R1,psth_R2,time,fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_cholinergic'));
    
    % Filtered HDB neurons
    avg_psth_cholinergic(selected_hdbcells, 'cue', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_hdb'));
    avg_psth_cholinergic(selected_hdbcells, 'reward', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_hdb'));
    avg_psth_cholinergic(selected_hdbcells, 'punish', 'none', fullfile(resdir, 'Fig2', PSTHdir1, 'filtered_hdb'));
    
    % Only include lick trials
    % All cholinergic neurons
    PSTHdir2 = 'avg_PSTH_lick_nonadaptive';
    avg_psth_cholinergic(achcells, 'lickcue', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'all_cholinergic'));
    avg_psth_cholinergic(achcells, 'lickreward', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'all_cholinergic'));
    avg_psth_cholinergic(achcells, 'lickpunish', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'all_cholinergic'));
    
    % All HDB neurons
    avg_psth_cholinergic(hdbcells, 'lickcue', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'all_hdb'));
    avg_psth_cholinergic(hdbcells, 'lickreward', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'all_hdb'));
    avg_psth_cholinergic(hdbcells, 'lickpunish', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'all_hdb'));
    
    % Filtered cholinergic neurons
    avg_psth_cholinergic(selected_achcells, 'lickcue', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'filtered_cholinergic'));
    avg_psth_cholinergic(selected_achcells, 'lickreward', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'filtered_cholinergic'));
    avg_psth_cholinergic(selected_achcells, 'lickpunish', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'filtered_cholinergic'));
    
    % Filtered HDB neurons
    avg_psth_cholinergic(selected_hdbcells, 'lickcue', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'filtered_hdb'));
    avg_psth_cholinergic(selected_hdbcells, 'lickreward', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'filtered_hdb'));
    avg_psth_cholinergic(selected_hdbcells, 'lickpunish', 'none', fullfile(resdir, 'Fig2', PSTHdir2, 'filtered_hdb'));
    
    % Omissions
    % Example
    PSTHdiro = 'omission_response';
    H = figure;
    cellid = 'HDB36_190426a_3.1';
    viewcell2b(cellid,'TriggerName','DeliverAllFeedback','SortEvent','TrialStart',...
        'sigma', 0.08,'eventtype','behav','ShowEvents',{{'StimulusOn'}},...
        'Partitions','#Omission','window',[-5 5])
    maximize_figure(H)
    
    cellidt = cellid;
    cellidt(cellidt=='.') = '_';
    fnm = fullfile(resdir,'Fig2',PSTHdiro,[cellidt '_O.jpg']);   % save
    saveas(H,fnm)
    close(H)
    
    % Average
    avg_psth_cholinergic(achcells, 'omission', 'none', fullfile(resdir, 'Fig2', PSTHdiro, 'all_cholinergic'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIGURE 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fig3spec
    
    % Reinforcement learning model
    modelresdir_allcells = fullfile(resdir, 'Fig3', 'all_cholinergic_bargraph_SD'); % results directory
    modelresdir_hdbcells = fullfile(resdir, 'Fig3', 'hdb');
    modelresdir_selected_achcells = fullfile(resdir, 'Fig3', 'selected_achcells');
    
    [a11, a21, S1, err1, c_a11, c_a21, c_S1, c_err1, model_values1] = cholinergic_RFmodel_main(achcells,...   % all neurons
        [-3 3],'baselinewin',[-3 -1],'responsewin_cue',[0 0.5],'responsewin_reinforcement',[0 0.2], 'resdir', modelresdir_allcells);
    %
        [a12, a22, S2, err2, c_a12, c_a22, c_S2, c_err2, model_values2] = cholinergic_RFmodel_main(hdbcells,... % hdb neurons only
            [-3 3],'baselinewin',[-3 -1],'responsewin_cue',[0 0.5],'responsewin_reinforcement',[0 0.2], 'resdir', modelresdir_hdbcells);
    
    [a13, a23, S3, err3, c_a13, c_a23, c_S3, c_err3, model_values3] = cholinergic_RFmodel_main(selected_achcells,... % hdb neurons only
        [-3 3],'baselinewin',[-3 -1],'responsewin_cue',[0 0.5],'responsewin_reinforcement',[0 0.2], 'resdir', modelresdir_selected_achcells);
    
    [H1, p1] = boxstat(err1,c_err1,'normal','reversed',0.05,'paired'); % boxplots to compare model fit
    [H2, p2] = boxstat(err2,c_err2,'normal','reversed',0.05,'paired');
    [H3, p3] = boxstat(a11,a21,'eta1','eta2',0.05,'paired'); % boxplot to compare model parameters for reward and punishment
    
    set(H1, 'Renderer', 'painters');
    set(H2, 'Renderer', 'painters');
    set(H3, 'Renderer', 'painters');

    fnm = fullfile(resdir, 'Fig3', 'model_stat_all_cholinergic.fig'); % save and close figures
    fnm2 = fullfile(resdir, 'Fig3', 'model_stat_hdb.fig');
    fnm3 = fullfile(resdir, 'Fig3', 'compare_params_all_cholinergic.fig');
    saveas(H1,fnm);
    close(H1)
    saveas(H2,fnm2);
    close(H2)
    saveas(H3,fnm3);
    close(H3)
    
    H1 = figure;
    boxplot(c_err1-err1); % boxplots to compare model fit
    H2 = figure;
    boxplot(c_err2-err2);
    fnm = fullfile(resdir, 'Fig3', 'model_stat_diff_all_cholinergic.fig'); % save and close figures
    fnm2 = fullfile(resdir, 'Fig3', 'model_stat_diff_hdb.fig');
    saveas(H1,fnm);
    close(H1)
    saveas(H2,fnm2);
    close(H2)
    
    % Simulate spiking data
    spikingsim(achcells,model_values1)
    
    % Correlation of model parameters with behavioral discimination
    resdir_lrc_allcells = fullfile(resdir,'Fig3','modellickratecorr','all_cholinergic'); % results directory
    modellickratecorr(achcells,a11,a21,resdir_lrc_allcells);
    %     resdir_lrc_hdbcells = fullfile(resdir,'Fig3','modellickratecorr','hdb'); % results directory
    %     modellickratecorr(hdbcells,a12,a22,resdir_lrc_hdbcells);
    
    resdir_lrc_selected_achcells = fullfile(resdir,'Fig3','modellickratecorr','selected_achcells'); % results directory
    modellickratecorr(selected_achcells,a13,a23,resdir_lrc_selected_achcells);
    
    % PSTH on simulated data
    PSTHdir1 = 'avg_PSTH_nonadaptive_sim_v4';
    avg_psth_cholinergic(achcells, 'cue', 'none', fullfile(resdir, 'Fig3', PSTHdir1, 'all_cholinergic'),'virtual')
    avg_psth_cholinergic(achcells, 'reward', 'none', fullfile(resdir, 'Fig3', PSTHdir1, 'all_cholinergic'),'virtual')
    avg_psth_cholinergic(achcells, 'punish', 'none', fullfile(resdir, 'Fig3', PSTHdir1, 'all_cholinergic'),'virtual')
    
    PSTHdir2 = 'avg_PSTH_nonadaptive_sim_v4_selected';
    avg_psth_cholinergic(selected_achcells, 'cue', 'none', fullfile(resdir, 'Fig3', PSTHdir2, 'filtered_cholinergic'),'virtual')
    avg_psth_cholinergic(selected_achcells, 'reward', 'none', fullfile(resdir, 'Fig3', PSTHdir2, 'filtered_cholinergic'),'virtual')
    avg_psth_cholinergic(selected_achcells, 'punish', 'none', fullfile(resdir, 'Fig3', PSTHdir2, 'filtered_cholinergic'),'virtual')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIGURE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if cue_rt    % Show if cue response is earlier than RT - FIG4.A-B
    resdir_cue_rt = fullfile(resdir, 'Fig4', 'cue_rt'); % results directory
    compare_cue_RT({'HDB18_170813a_4.2'}, resdir_cue_rt) % example
    compare_cue_RT({'HDB36_190510a_3.2'}, resdir_cue_rt) % example
    compare_cue_RT(achcells, resdir_cue_rt) % average of all cholinergic neurons
end

if cue_lick_corr % Lick response partitioned by cue response (high vs. low firing after cues, median split)
    resdir_corr = fullfile(resdir, 'Fig4', 'cue_lick_correlation'); % results directory
    cue_lick_correlation(achcells, resdir_corr);
end

if islickPSTH % average PSTHs partitioned by licking behavior (lick vs. no lick) - FIG4.C-E
    PSTHdir1 = 'avg_PSTH_nonadaptive';     % All cholinergic neurons
    avg_psth_cholinergic(achcells, 'cue', 'lick', fullfile(resdir, 'Fig4', PSTHdir1, 'all_cholinergic_islick'));
    avg_psth_cholinergic(achcells, 'reward', 'lick', fullfile(resdir, 'Fig4', PSTHdir1, 'all_cholinergic_islick'));
    avg_psth_cholinergic(achcells, 'punish', 'lick', fullfile(resdir, 'Fig4', PSTHdir1, 'all_cholinergic_islick'));
end

if rt_lick_corr % RT predicts cholinergic response and vice versa
    resdir_rtcurves1 = fullfile(resdir, 'Fig4', 'rtcurves_doubly_adaptive', 'all_cholinergic'); % results directory
    resdir_rtcurves2 = fullfile(resdir, 'Fig4', 'rtcurves_doubly_adaptive', 'hdb'); % results directory
    rtcurves_cholinergic(achcells,'align', 'tone', 'normalization', 'zscore', 'doraster', false, 'issave', true, 'resdir', resdir_rtcurves1);
    rtcurves_cholinergic(hdbcells,'align', 'tone', 'normalization', 'zscore', 'doraster', false, 'issave', true, 'resdir', resdir_rtcurves2);
    rtcurves_stats(fullfile(resdir, 'Fig4'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXTRA ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cholinergic PSTHs sorted by location
if psth_by_area
    areas = {'HDB' 'VP' 'MS' 'AcbC'};
    resdir = fullfile(getpref('cellbase', 'datapath'), '_adaptive_PSTH');
    for i = 1:length(areas)
        area = areas{i}; % current area for filter location
        cholinergic_PSTH(area, resdir)
    end
end

% CCG with other noncholinergic neurons (cholinergic-noncholinergic pairs)
if ccg_pairs
    ach_nonach_ccg(achcells, ccg_resdir,true);
end

%   Correlation of cue and reward response
if cue_reward_corr
    resdir_corr = fullfile(getpref('cellbase', 'datapath'), '_cue_reward_correlations');
    for k = 1:length(achcells)
        cCell = achcells{k}; % current cell
        response_correlation(cCell, 'cue', 'reward', resdir_corr)
    end
    
    resdir_c_r_corr1 = fullfile(resdir, 'Fig4', 'cue_reward_corr', 'all_cholinergic'); % results directory
    resdir_c_r_corr2 = fullfile(resdir, 'Fig4', 'cue_reward_corr', 'hdb'); % results directory
    cue_reward_correlation_cholinergic(achcells, resdir_c_r_corr1);
    cue_reward_correlation_cholinergic(hdbcells, resdir_c_r_corr2);
end

% PSTH conditioned on large/small cue response or previous outcome (reward/punishment)
if more_PSTHs
    PSTHdir1 = 'avg_PSTH_cueconditioned';
    avg_psth_cholinergic(achcells, 'cue', 'cue', fullfile(resdir, 'Fig4', PSTHdir1, 'all_cholinergic'))
    avg_psth_cholinergic(achcells, 'reward', 'cue', fullfile(resdir, 'Fig4', PSTHdir1, 'all_cholinergic'))
    avg_psth_cholinergic(achcells, 'punish', 'cue', fullfile(resdir, 'Fig4', PSTHdir1, 'all_cholinergic'))
    
    avg_psth_cholinergic(selected_achcells, 'cue', 'cue', fullfile(resdir, 'Fig4', PSTHdir1, 'filtered_cholinergic'))
    avg_psth_cholinergic(selected_achcells, 'reward', 'cue', fullfile(resdir, 'Fig4', PSTHdir1, 'filtered_cholinergic'))
    avg_psth_cholinergic(selected_achcells, 'punish', 'cue', fullfile(resdir, 'Fig4', PSTHdir1, 'filtered_cholinergic'))
    
    avg_psth_cholinergic(hdbcells, 'cue', 'cue', fullfile(resdir, 'Fig4', PSTHdir1, 'all_hdb'))
    avg_psth_cholinergic(hdbcells, 'reward', 'cue', fullfile(resdir, 'Fig4', PSTHdir1, 'all_hdb'))
    avg_psth_cholinergic(hdbcells, 'punish', 'cue', fullfile(resdir, 'Fig4', PSTHdir1, 'all_hdb'))
    
    history_variables(achcells)
    PSTHdir1 = 'avg_PSTH_historyconditioned';
    avg_psth_cholinergic(achcells, 'cue', 'previous_outcome', fullfile(resdir, 'Fig4', PSTHdir1, 'all_cholinergic'))
    avg_psth_cholinergic(achcells, 'reward', 'previous_outcome', fullfile(resdir, 'Fig4', PSTHdir1, 'all_cholinergic'))
    avg_psth_cholinergic(achcells, 'punish', 'previous_outcome', fullfile(resdir, 'Fig4', PSTHdir1, 'all_cholinergic'))
    
    avg_psth_cholinergic(hdbcells, 'cue', 'previous_outcome', fullfile(resdir, 'Fig4', PSTHdir1, 'all_hdb'))
    avg_psth_cholinergic(hdbcells, 'reward', 'previous_outcome', fullfile(resdir, 'Fig4', PSTHdir1, 'all_hdb'))
    avg_psth_cholinergic(hdbcells, 'punish', 'previous_outcome', fullfile(resdir, 'Fig4', PSTHdir1, 'all_hdb'))
end