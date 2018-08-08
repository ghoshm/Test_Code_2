%% ASD Cluster Frequency Figure

%% Load Data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'states');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'i_experiment_reps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'i_group_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'lb_merge');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'experiment_reps'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'fps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'cmap');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'night_color');

% Settings 
er = 1; % set experiment of interest
set_token = find(experiment_reps == er,1,'first'); % used for each experiments sets settings
time_bins = fps{set_token}*60*15; % set smoothing window 
m = [2 3]; % modules of interest (including inactive and active!).  

% Smooth Data 
for c = 1:length(m) % for each module of interest 
        clear data;
        
        data = states{er,1}(:,...
            lb_merge{er,1}(5):lb_merge{er,1}(7)) == m(c);
        
        tic 
        for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each fish
            smoothed_clusters{er,c}(g,:) = ...
                smooth(nanmean(data(i_group_tags(i_experiment_reps == er) == g,:)),time_bins); 
        end
        toc 
                
end

%% Cluster Frequencies Figure
% Run one er @ a time to avoid low level graphics errors

figure; % make a figure for this experiment
for c = 1:length(m) % for each module of interest 
    subplot(1,2,c); hold on;
    set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
    for g = size(smoothed_clusters{er,1},1):-1:1 % For each group
        legend_lines(1,g) = plot(lb_merge{er,1}(5):lb_merge{er,1}(7),...
            smoothed_clusters{er,c}(g,:),'color',cmap{set_token}(g,:),'linewidth',5);
        num(g) = sum(i_group_tags(i_experiment_reps == er) == g);
    end
    
    % Find the top & bottom
    scrap(1,1) = 0.5;
    scrap(2,1) = 0;
    
    % Night patches
    r = rectangle('Position',[lb_merge{er,1}(6) 0 (lb_merge{er,1}(7) - lb_merge{er,1}(6)) ... 
        max(scrap(:)) + max(scrap(:))*0.05],...
        'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
    uistack(r,'bottom'); % Send to back
    
    % Axis etc
    axis([(lb_merge{er,1}(5) + (1 + time_bins))...
        (lb_merge{er,1}(7) - time_bins) ...
        scrap(2,1) scrap(1,1)]);
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32);
    xlabel('Time (Hours)','Fontsize',32);
    set(gca,'XTick',[lb_merge{er,1}(6)]); 
    set(gca,'XTickLabels',{'23:00'}); 
    ylabel('Mean Frequency','Fontsize',32);
    
    % legend
    if c == length(m)
        legend(legend_lines,{horzcat('CHD8^{+/+},',' n = ',num2str(num(1))),...
            horzcat('CHD8^{-/+},',' n = ',num2str(num(2))),...
            horzcat('CHD8^{-/-},',' n = ',num2str(num(3)))});
        legend('Boxoff');
    end
    
end 

