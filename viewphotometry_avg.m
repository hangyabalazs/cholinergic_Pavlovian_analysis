function viewphotometry_avg(cellids, resdir, varargin)
%VIEWPHOTOMETRY_AVG   Average PSTH of fiber photometry data.
%   VIEWPHOTOMETRY_AVG(CELLIDS, RESDIR, VARARGIN) creates an average PSTH
%   of fiber photometry data recorded from sessions listed in CELLIDS.
%   PETHs are saved to RESDIR.

%   See also VIEWPHOTOMETRY.M, PERIEVENT.M

%   Balint Kiraly, Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu

% Default arguments
default_args={...
    'window',               [-3 3];...
    'dt',                   0.01;...
    'sigma',                10;...
    'isadaptive'            false;...
    'FigureNum',            1;...
    'Signal'                'dff' %s465, s405, dff
    'TriggerEvent',         'DeliverFeedback';...
    'SortEvent',            'TrialStart';...
    'ShowEvents',           {{'StimulusOn'}};...
    'ShowEventsColors',     {{[0 0.8 0] [0.8 0.8 0] [0 0.8 0.8]}};...
    'Num2Plot',             'all';...
    'PlotDashedEvent',      '';...
    'PlotDashedCondition',  'min';...
    'PSTHPlot',             1;...
    'PSTHlinewidth',        1.5;...
    'DashedLineStyle',      ':';...
    'LastEvents',           '';...
    'Partitions',           'all';...
    'PrintCellID',          'on';...
    'PrintCellIDPos',       'bottom-right';...
    'BurstPSTH'             'off';......
    };
[g,error] = parse_args(default_args,varargin{:});

% Input argument check
if nargin < 1
    error('Not enough input arguments. Please provide, at least, an animal and a session.');
end

isnorm = 0;

if ~isfolder(resdir)
    mkdir(resdir)
end

if strcmp(g.Partitions, '#AllReward') % baseline windows used for normalization
    bwin = [-3 -2];
    twin = 1;
    color1 = [0.22 0.78 0.62];
    color2 = [0 0.4 0.4];
elseif strcmp(g.Partitions, '#Reward') % baseline windows used for normalization
    bwin = [-3 -2];
    twin = 0.7;
    color1 = [0.22 0.78 0.62];
elseif strcmp(g.Partitions, '#Punishment')
    bwin = [-3 -2];
    twin = 1;
    color1 = [0.2 0 0];
    color2 = [1 0.6 0.8];
elseif strcmp(g.Partitions, '#Omission')
    bwin = [-3 -2];
    twin = 1;
    color1 = [0 0 0.2];
    color2 = [0 0 1];
elseif strcmp(g.Partitions, '#TrialType')
    bwin = [-1 0];
    twin = 0.5;
    color1 = [0.2 0.6 0.6];
    color2 = [0.4 0 0.2];
end

[psths1, psths1_st, psths2, psths2_st] = deal(nan(length(cellids),72289));

[mean1, mean2, mean1_st, mean2_st, maxval1, maxval2] = deal(nan(1,size(cellids,2)));


for c = 1:length(cellids) % Loop through sessions
    current_cellid = cellids{c};
    animalID = current_cellid{1};
    sessionID = current_cellid{2};
    channel = current_cellid{3};
    
    if strcmp(channel, 'Ch1')
        g.Signal = 'dff_D';
    elseif strcmp(channel, 'Ch2')
        g.Signal = 'dff_A';
    end
    % Load data and FiberEvents
    path = ['D:\HDB_cholinergic_reinf_learning\Fiber_photometry' filesep animalID filesep sessionID filesep];
    TE = load([path 'FiberEvents.mat']);
    DATA = load([path filesep 'proF.mat']);
    TE.Blocknum(TE.Hit~=1)=NaN;
    % Extracting valid trials
    [COMPTRIALS, TAGS] = partition_trials(TE,g.Partitions);
    vinx = cellfun(@(s)(~isempty(s)),COMPTRIALS);
    COMPTRIALS = COMPTRIALS(vinx);
    
    TAGS = TAGS(vinx);
    trigev = TE.(g.TriggerEvent);
    if ~iscell(trigev)
        valid_trials = find(~isnan(trigev));
    else
        valid_trials = find(cellfun(@(s)~isempty(s),trigev));
    end
    
    % Creating time vector
    time = g.window(1):(1/DATA.sr):g.window(2); % time vector
    b_inx = find((time>=bwin(1))&time<bwin(2)); % baseline window indices
    t_inx = find((time>=0)&(time<twin)); % test window indices
        
    % Sort trials
    NUMevents = length(g.SortEvent);
    if iscellstr(g.SortEvent)
        sort_var = nan(NUMevents,NUMtrials);
        for iS = 1:NUMevents
            sort_var(iS,:) = TE.(g.SortEvent{iS}) - TE.(g.TriggerEvent);
        end
        sort_var = min(sort_var);
    elseif ~isempty(g.SortEvent)
        if ~iscell(TE.(g.TriggerEvent))
            sort_var = TE.(g.SortEvent) - TE.(g.TriggerEvent);
        else
            gte = nan(1,NUMtrials);
            inx = ~cellfun(@isempty,TE.(g.TriggerEvent));
            gte(inx) = cell2mat(cellfun(@(s)s(1),TE.(g.TriggerEvent)(inx),'UniformOutput',false));
            sort_var = TE.(g.SortEvent) - gte;
        end
    else
        sort_var = NaN;
    end
    
    [mylabels, mycolors, mycolors2,mylinestyle] = makeColorsLabels(@defineLabelsColors_Balazs,TAGS);
    
    
    % Peri-event matrix for the given contingencies
    for iPAR = 1:length(TAGS)
        valT = intersect(valid_trials,COMPTRIALS{1,iPAR});
        sortT = TE.(g.TriggerEvent)(valT);
        vdisc = [];
        if strcmp(g.TriggerEvent,'TrialStart')
            [~,vdisc] = min(abs((sortT) - (DATA.tss-DATA.tss(1))));
        else
            startT = TE.TrialStart(valT);
            %vdisc = round((startT + sortT) * DATA.sr);
            for i= 1:length(startT)
                [~,vdisc(i)] = min(abs((startT(i) + sortT(i)) - (DATA.tss-DATA.tss(1))));
            end
        end
        
        nT = [1:1:length(vdisc)];
        
        [fibmean,fibcv,spmat] = perievent(vdisc,DATA.(g.Signal),DATA.sr,g.window,isnorm,g.dt,g.sigma); % peri-event
        
        
        [~,inx] = sort(sort_var(valT));
        
        if max(inx)>size(spmat,1)
            delinx = find(inx>size(spmat,1))
            inx(delinx)=[];
        end
        
        % Remove big artifacts
        spmat_avg = mean(spmat, 1);
        spmat_std = std(spmat);
        boundary_up = spmat_avg+2*spmat_std;
        boundary_low = spmat_avg-2*spmat_std;
        
        for j = 1:size(inx,2)
            all_points = length(boundary_up);
            ten_pct = all_points/10;
            deviation_up=find(spmat(inx(j),:)>boundary_up);
            deviation_low = find(spmat(inx(j),:)<boundary_low);
            if (length(deviation_up)+length(deviation_low))>ten_pct
                inx(j)=nan;
            end
        end
        inx(isnan(inx)) = []; % Only keep low noise trials
        
        if iPAR == 1
            psths1(c,:) = nanmean(spmat(inx,:),1);
            mn = nanmean(psths1(c,b_inx));
            sd = nanstd(psths1(c,b_inx));
            psths1_st(c,:) = (psths1(c,:)-mn)/sd;
        elseif iPAR == 2
            psths2(c,:) = nanmean(spmat(inx,:),1);
            psths2_st(c,:) = (psths2(c,:)-mn)/sd; % use the same mean and sd for normalization
        end
    end
end

for m = 1:size(psths1,1)
    mean1(m) = nanmean(psths1(m,t_inx));
    mean2(m) = nanmean(psths2(m,t_inx));
    
    mean1_st(m) = nanmean(psths1_st(m,t_inx));
    mean2_st(m) = nanmean(psths2_st(m,t_inx));
    
    maxval1(m)= max(psths1(m,t_inx));
    maxval2(m)= max(psths2(m,t_inx));
end

tag = g.Partitions(2:end);

% Boxstat without normalization
boxstat2(mean1, mean2, 'T1', 'T2', 0.005, 'paired')
err=se_of_median(mean1-mean2);
fnm = fullfile(resdir,['Boxstat' tag '.fig']);
saveas(gcf,fnm)
set(gcf,'renderer', 'painters')
fnm = fullfile(resdir,['Boxstat' tag '.eps']);
saveas(gcf,fnm)
close(gcf)

% Boxstat with normalized values
boxstat2(mean1_st, mean2_st, 'T1', 'T2', 0.005, 'paired')
err_st=se_of_median(mean1_st-mean2_st);
fnm = fullfile(resdir,['Boxstat_normalized_' tag '.fig']);
saveas(gcf,fnm)
set(gcf,'renderer', 'painters')
fnm = fullfile(resdir,['Boxstat_normalized_' tag '.eps']);
saveas(gcf,fnm)
% close(gcf)

% Bar graph of mean difference w/o normalization
figure;
bar([1],[median(mean1-mean2)])
hold on
errorbar([1],[median(mean1-mean2)],err,err);
fnm = fullfile(resdir,['Barplot_' tag '.fig']);
saveas(gcf,fnm)
set(gcf,'renderer', 'painters')
fnm = fullfile(resdir,['Barplot_' tag '.eps']);
saveas(gcf,fnm)
close(gcf)

% Bar graph of mean difference with normalization
figure;
bar([1],[median(mean1_st-mean2_st)])
hold on
errorbar([1],[median(mean1_st-mean2_st)],err_st,err_st);
fnm = fullfile(resdir,['Barplot_norm_' tag '.fig']);
saveas(gcf,fnm)
set(gcf,'renderer', 'painters')
fnm = fullfile(resdir,['Barplot_norm_' tag '.eps']);
saveas(gcf,fnm)
close(gcf)

figure;
boxplot(mean1_st-mean2_st);
hold on;
y = mean1_st - mean2_st;
x = ones(length(y),1);
scatter(x(:),y(:),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.05);
I = gcf;
fnm = fullfile(resdir,['Boxplot_normalized_' tag '.fig']);
saveas(I,fnm)
set(I,'renderer', 'painters')
fnm = fullfile(resdir,['Boxplot_normalized_' tag '.eps']);
saveas(I,fnm)

figure;
errorshade(time, nanmean(psths1_st,1), nanstd(psths1_st,1)/sqrt(size(psths1_st,1)),...
    'LineColor',color1,'ShadeColor',color1)
hold on
errorshade(time, nanmean(psths2_st,1), nanstd(psths2_st,1)/sqrt(size(psths2_st,1)),...
    'LineColor',color2,'ShadeColor',color2)
hold off
fnm = fullfile(resdir,['avg_PETH_' tag '.fig']);
saveas(gcf,fnm)
set(gcf,'renderer', 'painters')
fnm = fullfile(resdir,['avg_PETH_' tag '.eps']);
saveas(gcf,fnm)
fnm = fullfile(resdir,['avg_PETH_' tag '.jpeg']);
saveas(gcf,fnm)
