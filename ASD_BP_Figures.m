%% ASD Bout Proportion Figures 

% Load Data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'bout_proportions');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'cmap');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'i_experiment_reps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'i_group_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731.mat', 'night_color');

%% Figure 
er = 2; % set experiment of interest 
set_token = 1; % setttings 
s = 2; % set state of interest 
m = [1 5]; % set modules of interest 

for c = 1:length(m) 
    subplot(1,2,c);
    set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
    clear data; 
    data = squeeze(bout_proportions{s,1}(i_experiment_reps == er,m(c),3:6)); 
    sep = 0; 
    for g = 1:max(i_group_tags(i_experiment_reps == er))
        spread_cols(g,:) = plotSpread(data(i_group_tags(i_experiment_reps == er) == g,:),...
            'XValues',[1:size(data,2)] + sep,'spreadWidth',.3,'distributionColors',...
            cmap{set_token}(g,:),'ShowMM',2);
        set(findall(gca,'type','line'),'markersize',15); % change
        spread_cols{g,2}.LineWidth = 3; spread_cols{g,2}.Color = [0 0 0]; % Change marker properties
        data_m(g,:) = nanmean(data(i_group_tags(i_experiment_reps == er) == g,:));
        sep = sep + .3;  
        num(g) = sum(i_group_tags(i_experiment_reps == er) == g); 
    end
    
    sep = 1:0.3:1.6; % hard coded 
    for t = 1:size(data,2) % for each time window 
        plot(sep+(t-1), data_m(:,t),'k','linewidth',2);
    end 
    
    a = 1; night_start = 2; 
    for n = 1:2 % For each night
        r(a) = rectangle('Position',[(night_start-0.2) 0 1 1],...
            'FaceColor',night_color{1},'Edgecolor',[1 1 1]);
        uistack(r(a),'bottom'); % Send to back
        a = a + 1; night_start = night_start + 2; % Add to counters
    end
    
    axis([0.8 4.8 0 .5]);
    set(gca,'XTick',[]);
    xlabel('Time (Days/Nights)','Fontsize',32);
    ylabel('Probability','Fontsize',32);

    if c == 2
        temp = get(gca,'Children');
        temp = [temp(16) ; temp(11) ; temp(6)]; % hard coded 
        legend(temp,{horzcat('BCKDK^{+/+},',' n = ',num2str(num(1))),...
            horzcat('BCKDK^{-/+},',' n = ',num2str(num(2))),...
            horzcat('BCKDK^{-/-},',' n = ',num2str(num(3)))});
        legend('boxoff');
    end
end 