%% Paper Feedback 

% Code for working on paper feedback 

%% Split-Matrix Figure 
        
% Load data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat', 'states');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat', 'fps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat', 'cmap_cluster');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat', 'lb_merge');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat', 'time_window'); 
time_bins = fps{1}*60*5; % set smoothing window (note assumes a constant frame rate across experiments)

% Sorting by activity 
temp = states{1,1}; 
temp(temp <= 5) = NaN; % hard coded
temp = temp(:,lb_merge{1,1}(time_window{1}(1)):lb_merge{1,1}(time_window{1}(2)+1)); 
temp(isnan(temp)==0) = 1; 
[~,O] = sort(nansum(temp,2)); 
states{1,1} = states{1,1}(flip(O),:);

% Active 
figure;
subplot(1,2,1); hold on; set(gca,'FontName','Calibri'); set(gca,'Fontsize',32);
title(horzcat('Active Modules'));
temp = states{1,1}; 
temp(temp <= 5) = NaN; % hard coded
temp = temp(:,lb_merge{1,1}(time_window{1}(1)):lb_merge{1,1}(time_window{1}(2)+1)); 
ax = imagesc(temp,...
    'AlphaData',isnan(temp)==0); 

colormap(gca,cmap_cluster{1,1}); % active colormap
box off; set(gca,'Layer','top'); 
set(ax,'CDataMapping','direct');
set(gca,'YDir','Reverse'); 
c = colorbar; c.Label.String = 'Module';
axis([(1 + 1 + time_bins)...
    ((size(temp,2) + 1) - time_bins) ...
    1 size(states{1,1},1)]);
xlabel('Time (Days/Nights)','Fontsize',32);
set(gca, 'XTick', []);
ylabel('Fish ID','Fontsize',32);

clear ax c x

% Inactive  
subplot(1,2,2); hold on; set(gca,'FontName','Calibri'); set(gca,'Fontsize',32);
title(horzcat('Inactive Modules'));
temp = states{1,1}; 
temp(temp > 5) = NaN; % hard coded
temp = temp(:,lb_merge{1,1}(time_window{1}(1)):lb_merge{1,1}(time_window{1}(2)+1)); 
ax = imagesc(temp,...
    'AlphaData',isnan(temp)==0); 

colormap(gca,cmap_cluster{2,1}); % inactive colormap
box off; set(gca,'Layer','top'); 
set(ax,'CDataMapping','direct');
set(gca,'YDir','Reverse'); 
c = colorbar; c.Label.String = 'Module';
axis([(1 + 1 + time_bins)...
    ((size(temp,2) + 1) - time_bins) ...
    1 size(states{1,1},1)]);
xlabel('Time (Days/Nights)','Fontsize',32);
set(gca, 'XTick', []);
ylabel('Fish ID','Fontsize',32);

clear ax c x

%% Real vs shuffled Library Size 

% Load data
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'grammar_size');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'i_experiment_reps'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'i_experiment_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'i_group_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'experiment_reps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'cmap');

%% Size Figure  
er = 1; % for the WT fish 

figure;
subplot(1,2,1);
set_token = find(experiment_reps == er,1,'first'); % settings
counter = 1; % counts groups for plots
hold on; set(gca,'FontName','Calibri'); clear scrap;

for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
    clear data;
    data = [repmat(grammar_size(i_experiment_reps == er & i_group_tags == g,1),(size(grammar_size,2)-1),1) ...
        reshape(grammar_size(i_experiment_reps == er & i_group_tags == g,2:end),[],1)];
    plot([counter,counter+1],data,...
        'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
    errorbar([counter,counter+1],nanmean(data),nanstd(data),...
        'color',cmap{set_token}(g,:),'linewidth',3);
    counter = counter + 2; % add to counter
    
    scrap(1,g) = min(min(grammar_size(i_experiment_reps == er & i_group_tags == g,:)));
    scrap(2,g) = max(max(grammar_size(i_experiment_reps == er & i_group_tags == g,:)));
end

box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
if er == 1 % for the WT Data
    set(gca,'XTick', [1 2]); % set X-ticks
    set(gca,'XTickLabels',{'Data','Shuffled'}); % X Labels
else % for the other experiments
    set(gca,'XTickLabels',geno_list{set_token}.colheaders); % X Labels
end
ylabel('Number of Motifs','Fontsize',32); % Y Labels
axis([.5 (max(i_group_tags(i_experiment_reps == er))*2)+1 ...
    (min(scrap(1,:)) - (min(scrap(1,:))*0.05)) (max(scrap(2,:)) + (max(scrap(2,:))*0.05))]);

% Insert 
axes('Position',[0.4 0.5 0.1 0.4]); hold on; 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',16); set(gca,'FontName','Calibri'); % Set Font
spread_cols = plotSpread(grammar_size(i_experiment_reps == er,1) - ...
    nanmean(grammar_size(i_experiment_reps == er,2:end),2),...
    'distributionColors',cmap{1}(1,:),'showMM',2);
spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change marker properties
set(findall(gca,'type','line'),'markersize',15); % change
ylabel('? Number of Motifs'); 
set(gca,'XTick',[]); 

%% Module Dists   

% Load Data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'uniqueSeqs');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'seq_lengths_pd');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'numComp'); 

% Construct Grammar Matricies  
for tc = 1:size(uniqueSeqs,2)
    grammar_mat{1,tc} = nan(size(uniqueSeqs{1,tc},1),size(seq_lengths_pd,2)+1,'single'); % sequences x max length
    
    for s = 1:size(uniqueSeqs{1,tc},1) % for each sequence
        grammar_mat{1,tc}(s,1:size(uniqueSeqs{1,tc}{s,1},2)) = uniqueSeqs{1,tc}{s,1}; % fill in sequence
    end
    
end

% Length Distributions
for tc = 1:size(uniqueSeqs,2)
    subplot(4,3,tc); hold on;
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',24); set(gca,'FontName','Calibri'); % Set Font
    plot([(numComp(2)+.5) (numComp(2)+.5)],[0 .25],'--k','linewidth',1.5);
    
    scrap = [sum(isnan(grammar_mat{1,tc})==0,2) grammar_mat{1,tc}]; % lengths & grammar
    for i = unique(scrap(:,1))' % for each motif length
        if length(find(scrap(:,1) == i)) < 3 % if there are less than 3 motifs of this length
            scrap(scrap(:,1) == i,:) = []; % remove them
        end
    end
    
    scrap_cmap = lbmap(length(unique(scrap(:,1))),'BlueGray'); % colormap
    counter = 1; % start counter
    for i = unique(scrap(:,1))' % for each motif length
        toy = scrap(scrap(:,1) == i,2:end); % take motifs of this length
        toy = toy(:); % vectorise
        toy(isnan(toy)) = []; % remove nan values
        legend_cols(counter) = plot(histcounts(toy,'BinMethod','integers','normalization','probability'),...
            'color',scrap_cmap(counter,:),'linewidth',3);
        counter = counter + 1;
    end
    
%     legend(legend_cols,string(unique(scrap(:,1))),...
%         'Location','northeast');
%     legend('boxoff');
    axis tight;
    xlabel('Module','Fontsize',24);
    ylabel('Probability','Fontsize',24);
end