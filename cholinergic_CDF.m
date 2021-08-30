function cholinergic_CDF(T1_vals, T2_vals, Sim_vals, resdir)
%CHOLINERGIC_CDF(T1_VALS, T2_VALS, SIM_VALS, RESDIR) Calculates and plots cumulative density function.
%   CHOLINERGIC_CDF(T1_VALS, T2_VALS, SIM_VALS, RESDiR) calculates and plots
%   cumulative density function and scatter plot for recorded neural
%   responses and simulated neural responses to reward predicting cue,
%   punishment predicting cue, expected and surprising reward and
%   punishment respectively. The results are saved to RESDIR.
%
%   See also HISTC, STAIRS

%   Hegedus Panna
%   16-Nov-2020
%   panna.hegedus@koki.hu

%   Code review: BH 11/20/20

% Convert simulated values into a better format
[cue_sim, rew_sim, pun_sim] = deal(nan(2,size(T1_vals,1)));
for k = 1:length(Sim_vals) % reorganize simulated values
    for j = 1:2
        cue_sim(j,k) = Sim_vals{k}(1,j);
        rew_sim(j,k) = Sim_vals{k}(2,j);
        pun_sim(j,k) = Sim_vals{k}(3,j);
    end
end

% Plot CDF and scatter plot for events (real and simulated FR)
H1 = figure;
H2 = figure;
colors = {[0 0 1] [0 1 0] [1 0 0]};
sim_val2plot = {cue_sim rew_sim pun_sim};
for a = 1:length(sim_val2plot)
    color_choice = colors{a};
    value_choice = sim_val2plot{a};
    figure(H1)
    subplot(2,3,a) % CDF of reward cue, expected reward and surprising punishment response
    edges = 0:0.0001:max(T1_vals(:,a)')+0.01;   % bin edges
    dist_spike_avg = histc(T1_vals(:,a)',edges);   % histogram
    dist_spike_avg = [0 dist_spike_avg(1:end-1)];   % values corresponding to the edges
    dist_spike_avg = dist_spike_avg / sum(dist_spike_avg);   % normalize
    stairs(edges,cumsum(dist_spike_avg),'Color',color_choice,'LineWidth',2)
    axis tight
    set(gca,'box','off','FontSize',12,'TickDir','out')
    setmyplot_balazs
    set(gcf, 'Renderer', 'painters')
    
    hold on % plot simulated values
    dist_spike_avg_sim = histc(value_choice(1,:),edges);   % histogram
    dist_spike_avg_sim = [0 dist_spike_avg_sim(1:end-1)];   % values corresponding to the edges
    dist_spike_avg_sim = dist_spike_avg_sim / sum(dist_spike_avg_sim);   % normalize
    stairs(edges,cumsum(dist_spike_avg_sim),'Color',[128 128 128]/255,'LineWidth',2)
    hold off
    
    figure(H2)
    subplot(2,3,a) % scatter plot of reward cue, expected reward and surprising punishment response
    scatter(T1_vals(:,a),value_choice(1,:))
    axis equal
    line([0 150],[0 150])
    axis tight
    
    figure(H1)
    subplot(2,3,3+a) % CDF of punish cue, surprising reward and expected punishment response
    edges = 0:0.0001:max(T2_vals(:,a)')+0.01;   % bin edges
    dist_spike_avg2 = histc(T2_vals(:,a)',edges);   % histogram
    dist_spike_avg2 = [0 dist_spike_avg2(1:end-1)];   % values corresponding to the edges
    dist_spike_avg2 = dist_spike_avg2 / sum(dist_spike_avg2);   % normalize
    stairs(edges,cumsum(dist_spike_avg2),'Color',color_choice,'LineWidth',2)
    axis tight
    set(gca,'box','off','FontSize',12,'TickDir','out')
    setmyplot_balazs
    set(gcf, 'Renderer', 'painters')
    
    hold on % plot simulated values
    dist_spike_avg_sim2 = histc(value_choice(2,:),edges);   % histogram
    dist_spike_avg_sim2 = [0 dist_spike_avg_sim2(1:end-1)];   % values corresponding to the edges
    dist_spike_avg_sim2 = dist_spike_avg_sim2 / sum(dist_spike_avg_sim2);   % normalize
    stairs(edges,cumsum(dist_spike_avg_sim2),'Color',[128 128 128]/255,'LineWidth',2)
    hold off
    
    figure(H2)
    subplot(2,3,3+a) % scatter plot of punish cue, surprising reward and expected punishment response
    scatter(T2_vals(:,a),value_choice(2,:))
    axis equal
    line([0 150],[0 150])
    axis tight
end

% Save figures
saveas(H1, fullfile(resdir, 'sim_responses_CDF.eps')) % save plot
saveas(H1, fullfile(resdir, 'sim_responses_CDF.jpg')) % save
saveas(H1, fullfile(resdir, 'sim_responses_CDF.fig')) % save
saveas(H2, fullfile(resdir, 'sim_responses_scatterplot.eps')) % save plot
saveas(H2, fullfile(resdir, 'sim_responses_scatterplot.jpg')) % save
saveas(H2, fullfile(resdir, 'sim_responses_scatterplot.fig')) % save