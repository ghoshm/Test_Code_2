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
