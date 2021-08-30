function rtcurves_stats(sourcedir)


load(fullfile(sourcedir, 'allpsth_data.mat')); % load previously saved data

inx = find(time>0 & time<=0.5); % time indices

allspsth_c500=allspsth(:,:,inx);
allspsth_c500_q1=mean(allspsth(1,:,inx),3);
allspsth_c500_q2=mean(allspsth(2,:,inx),3);
allspsth_c500_q3=mean(allspsth(3,:,inx),3);
allspsth_c500_q4=mean(allspsth(4,:,inx),3);

mean_cueresp = [mean(allspsth_c500_q1) mean(allspsth_c500_q2) mean(allspsth_c500_q3) mean(allspsth_c500_q4)];
sem_cueresp = [std(allspsth_c500_q1)/5 std(allspsth_c500_q2)/5 std(allspsth_c500_q3)/5 std(allspsth_c500_q4)/5]



figure; % bar plot
bar([1 2 3 4], mean_cueresp);
hold on
errorbar([1 2 3 4], mean_cueresp, sem_cueresp)


% ANOVA
[p, table]=anova1([allspsth_c500_q1' allspsth_c500_q2' allspsth_c500_q3' allspsth_c500_q4'])