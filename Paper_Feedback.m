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

%% Analysis Example 

% Load Data
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\mRMR_Comparisons\Comps_250\Melatonin_D.mat');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Grammar_Results_Final.mat', 'gCount');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'grammar_mat');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'numComp');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'sleep_cells');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'idx_numComp_sorted');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'cmap_cluster_merge');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'bouts');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'i_experiment_reps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'i_group_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'cmap');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'geno_list');

% Choose dataset 
er = 3; 

% Choose Motif (Melatonin 44,768)
s = comps_v{er,1}(1,1);

% mean inactive module length (frames)
ibl = grpstats(sleep_cells(:,3),idx_numComp_sorted{2,1},'mean');
ibl(1) = []; % remove NaN's
ibl(end) = 250; % crop longest module to 10s (hard coded)

% Plot the Motif
figure; hold on; set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
b = 1; % baseline counter (plots from bottom to top)
seq = grammar_mat{1,1}(s,:); % find this motifs module sequence
seq(isnan(seq)) = []; % remove nan values
a = 1; % start a counter (frames)

for t = 1:length(seq) % for each module in the sequence
    if seq(t) <= numComp(1) % for the inactive modules
        plot([a (a+ibl(seq(t)))],[b b],...
            'color',cmap_cluster_merge(seq(t),:),'linewidth',15); % plot
        a = a + ibl(seq(t)); % add to time
    else % for the active modules
        plot(a:(a+length(nanmean(bouts{1,seq(t)-numComp(1)}))+1),...
            [b ((nanmean(bouts{1,seq(t)-numComp(1)})/28)+b) b],...
            'color',cmap_cluster_merge(seq(t),:),'linewidth',15); % plot
        a = a + length(nanmean(bouts{1,seq(t)-numComp(1)})) + 1; % add to time
    end
end

axis tight
axis off

% Count motif usage  
a = 1; 
for f = find(i_experiment_reps == er)'
    data(1,a) = gCount{f,1}(s,1);
    a = a + 1;
end 

% Figure 
figure; hold on; set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group 
    spread_cols = plotSpread(data(:,i_group_tags(i_experiment_reps == er) == g)',...
        'xyOri','flipped','spreadWidth',1/max(i_group_tags(i_experiment_reps == er)),...
        'distributionColors',...
        cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'XValues',...
        (1:size(data,1)) + (max(i_group_tags(i_experiment_reps == er))-1)/10 - (g-1)/10);
    set(findall(gca,'type','line'),'markersize',30); % change marker size 
    
    legend_cols(g,:) = errorbar(nanmean(data(:,i_group_tags(i_experiment_reps == er) == g),2),...
        ((1:size(data,1)) + (max(i_group_tags(i_experiment_reps == er))-1)/10 - (g-1)/10),...
        nanstd(data(:,i_group_tags(i_experiment_reps == er) == g)'),...
        'horizontal','o','linewidth',3,'color',cmap{set_token}(g,:),'capsize',9);
    
end


xlabel('Motif Counts','Fontsize',32);
set(gca,'YTick',[]);
axis([minmax(data) ylim])

%% Temporal Correlation (Seconds)

% Settings 
secs = 15; % hard coded, correlation lag window (seconds)
    % NOTE: Must look for longer than 20.04 seconds 
    
% Load data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat', 'states');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat', 'fps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat', 'lb_merge');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat', 'time_window'); 

temp = states{1,1}; % fish x frames
temp = temp(:,lb_merge{1,1}(time_window{1}(1)):lb_merge{1,1}(time_window{1}(2)+1)); % crop data  

% Pre-allocate 
scrap_corr = NaN(size(temp,1),(secs*fps{1}*2)+1,100,'single'); % fish x corr lags x pariwise comparisons (hard coded)

for f = 1 % for each fish
    tic
    scrap_corr_f = temp(f,:); % fish f's data
    counter = 1; % counts pairwise comparisons
    
    for a = 1:10 % for each module 
        for b = 1:10 % for each module  
            scrap_corr(f,:,counter) = ...
                xcorr(subplus(diff(scrap_corr_f == a)),subplus(diff(scrap_corr_f == b)),...
                secs*fps{1},'coeff');
            % 'coeff' — Normalizes the sequence so that the autocorrelations at zero lag equal 1:
            counter = counter + 1; % add to counter 
        end 
    end 
    toc
    
end 

%% Coding the transitions 
clear scrap_corr_trans
scrap_corr_trans = [ones(5,5) ; ones(5,5)*2 ]; 
scrap_corr_trans = [scrap_corr_trans scrap_corr_trans + 2]; 
scrap_corr_trans = scrap_corr_trans(:); 

%% Figure 
% NOTE: may be good to smooth the graphs
cols = NaN(length(scrap_corr_trans),3); 

for c = 1:4 
    cols(scrap_corr_trans == c,:) = lbmap(25,'BlueGray'); 
end 

clf;
for i = 1:4
subplot(1,4,i); axis([-5,5,0,0.07]); hold on; 
end
for i = 1:100
    subplot(1,4,scrap_corr_trans(i)); 
    plot(-secs:(1/fps{1}):secs,scrap_corr(1,:,i),'color',cols(i,:));
end

% Insert
axes('Position',[0.8 0.5 0.1 0.4]); hold on; 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',16); set(gca,'FontName','Calibri'); % Set Font
counter = 1; % start a counter 
for a = 1:10 % for each module
    for b = 1:10 % for each module 
        scatter(a,b,90,'markerfacecolor',cols(counter,:),...
            'markeredgecolor',cols(counter,:)); % scatter a marker
        counter = counter + 1; % add to counter 
    end
end

%% Temporal Correlation (Modules)

% Settings 
mods = 50; % corr lags 
er = 1; % experiments of interest 

% Load data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'threads');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat', 'time_window');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'i_experiment_reps');

% Pre-allocate 
scrap_corr = NaN(sum(i_experiment_reps == 1),(mods*2+1),100,'single'); % fish x corr lags x pariwise comparisons (hard coded)

for f = 1 % for each fish
    tic
    scrap_corr_f = threads{1,1,1}(threads{1,3,1} >= 3 & threads{1,3,1} <= 6); % fish f's data
    counter = 1; % counts pairwise comparisons
    
    for a = 1:10 % for each module 
        for b = 1:10 % for each module  
            scrap_corr(f,:,counter) = ...
                xcorr(scrap_corr_f == a,scrap_corr_f == b,...
                mods,'coeff');
            % 'coeff' — Normalizes the sequence so that the autocorrelations at zero lag equal 1:
            scrap_corr_labs(counter,1:2) = [a b];
            counter = counter + 1; % add to counter 
        end 
    end 
    toc
    
end 

%% Coding the transitions 
clear scrap_corr_trans
scrap_corr_trans = [ones(5,5) ; ones(5,5)*2 ]; 
scrap_corr_trans = [scrap_corr_trans scrap_corr_trans + 2]; 
scrap_corr_trans = scrap_corr_trans(:); 
scrap_corr_trans(scrap_corr_trans == 3) = NaN; 
scrap_corr_trans(scrap_corr_trans == 4) = 3; 

%% Figure 
test = squeeze(scrap_corr(1,:,:))';
test(test <= 0) = NaN;

% Single Panel 
figure;
data = test(isnan(scrap_corr_trans) == 0,:); % take transitions of this type
data_t = scrap_corr_labs(isnan(scrap_corr_trans) == 0,:);
[~,idx] = sort(data_t(:,1),'descend'); % sort by total correlation
data = data(idx,:); % sort by total correlation
data_t = data_t(idx,:); % sort by total correlation

subplot(1,2,1);
ax = imagesc(data_t,'AlphaData',isnan(data_t)==0); % plot data
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]);
set(ax,'CDataMapping','direct');

ax = subplot(1,2,2);
imagesc(data,'AlphaData',isnan(data)==0); % plot data
colormap(ax,'default')

% Multiple Panels 
clf;
for i = 1:3 % for each type of transition
    subplot(1,3,i); % subplot
    data = test(scrap_corr_trans == i,:); % take transitions of this type
    data_t = scrap_corr_labs(scrap_corr_trans == i,:);
    [~,idx] = sort(nansum(data,2),'descend'); % sort by total correlation
    data = data(idx,:); % sort by total correlation 
    data_t = data_t(idx,:); % sort by total correlation
    ax = imagesc(data_t,'AlphaData',isnan(data_t)==0); % plot data 
    colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]);
    set(ax,'CDataMapping','direct');
end

% Insert
axes('Position',[0.8 0.5 0.1 0.4]); hold on; 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',16); set(gca,'FontName','Calibri'); % Set Font
counter = 1; % start a counter 
for a = 1:10 % for each module
    for b = 1:10 % for each module 
        scatter(a,b,90,'markerfacecolor',cols(counter,:),...
            'markeredgecolor',cols(counter,:)); % scatter a marker
        counter = counter + 1; % add to counter 
    end
end

%% Compressibility vs Shuffle 
    % Based on Compression_values_test.m 
    
% Load data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Compression_Values_Test.mat');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'experiment_reps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'i_group_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'i_experiment_reps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'cmap');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'geno_list');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'i_experiment_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'nights');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'night_color');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'time_window');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'days_crop');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'days');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'nights_crop');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'nights');

%% Remove Spare edges 

for f = 1:size(threads,1) % for each fish 
    
    chunks{f,1}(end) = []; % remove the spare edge
    
end 
 
%% Reshape Chunks & Calculate Compression  
    % The compressibility of a sequence of uncompressed length l is given by the sum of the savings S
    % at each iteration divided by l (Gomez-Marin et al.,2016)
    
% settings
tw = 48; % hard coded maximum number of time windows
dn_hour = ones(1,tw)*2; % day and night hours
dn_hour([1:14 25:38]) = 1; % day (1) and night (2) hours

% allocate
compressibility = nan(size(threads,1),tw,size(threads,3),'single'); % fish x max time hour bins x tc

% fill data
for f = 1:size(threads,1) % for each fish
        
    for tc = 1:size(threads,3) % for real and shuffled data
        for h = 1:max(threads{f,3,1}) % for each hour
            try
                compressibility(f,h,tc) = nanmean(totSavings_cells{f,1}(tc,...
                    threads{f,3,1}(chunks{f,1}) == h)/step);
                % take a mean score every hour
            catch
                % catches hours with less cluster than chunks
                % note that these will remain NaN
            end
        end
    end
end

clear f h 

%% Calculate Difference from Shuffled Every Hour 

compressibility = compressibility(:,:,1) - ...
    nanmean(compressibility(:,:,2:end),3); 

%% %% Compressibility Every Hour - Figure
er = 1; 
set_token = find(experiment_reps == er,1,'first'); % settings
figure;
hold on; set(gca,'FontName','Calibri'); clear scrap;

for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
    clear data;
    data = compressibility(i_experiment_reps == er & i_group_tags == g,:);
    
    if er == 1 % for the WT experiments
        plot(data',...
            'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
        scrap(1,g) = min(data(:));
        scrap(2,g) = max(data(:));
        
    else 
        scrap(1,g) = min(nanmean(data) - nanstd(data));
        scrap(2,g) = max(nanmean(data) + nanstd(data)); 
    end
    
    legend_lines(g) = errorbar(nanmean(data),nanstd(data),...
        'color',cmap{set_token}(g,:),'linewidth',3);
    legend_cell{g} = horzcat(geno_list{set_token}.colheaders{g},', n = ',...
            num2str(size(data,1)));
end

y_lims = [(min(scrap(1,:)) + min(scrap(1,:))*0.05) ...
    (max(scrap(2,:)) + max(scrap(2,:))*0.05)]; % Add a bit of space either side
   
% Night Patches 
a = 1; night_start = 15; % hard coded counter
for n = 1:size(nights{set_token},2) % For each night
    r(a) = rectangle('Position',[(night_start) y_lims(1)...
        9 (y_lims(2)-y_lims(1))],...
        'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
    uistack(r(a),'bottom'); % Send to back
    a = a + 1; night_start = night_start + 24; % Add to counters
end
    
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
xlabel('Time (Hours)','Fontsize',32); % X Labels 
ylabel({'? Compressibility' ; '(per 500 modules)'},'Fontsize',32); % Y Labels
axis([1 (n*24) y_lims]); 
if er == 1 % for the WT data
   [~,icons,plots,~] = legend(legend_lines,legend_cell,'Location','northwest');
else
   [~,icons,plots,~] = legend(legend_cell,'Location','northeast');
end
legend('boxoff'); 
set(icons(1:g),'Fontsize',32) ; set(plots,'LineWidth',3);

clear et set_token g data scrap legend_lines legend_cell y_lims a night_start n r icons plots   

%% Compressibility Day vs Night
% settings
tw = 48; % hard coded maximum number of time windows
dn_hour = ones(1,tw)*2; % day and night hours
dn_hour([1:14 25:38]) = 1; % day (1) and night (2) hours

figure;
for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % settings
    subplot(2,2,er); counter = 1; % counts groups for plots
    hold on; set(gca,'FontName','Calibri'); clear scrap;
    
    for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
        clear data;
        % average day per fish
        data(:,1) = nanmean(compressibility(i_experiment_reps == er & i_group_tags == g,dn_hour == 1),2);
        % average night per night
        data(:,2) = nanmean(compressibility(i_experiment_reps == er & i_group_tags == g,dn_hour == 2),2);
        plot([counter,counter+1],data,...
            'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
        errorbar([counter,counter+1],nanmean(data),nanstd(data),...
            'color',cmap{set_token}(g,:),'linewidth',3);
        counter = counter + 2; % add to counter
        
        scrap(1,g) = min(data(:));
        scrap(2,g) = max(data(:));
    end
    
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
    if er == 1 % for the WT Data
        set(gca, 'XTick', [1 2]); % set X-ticks
        set(gca,'XTickLabels',{'Day','Night'}); % X Labels
    else % for the other experiments
        set(gca,'XTick',1.5:2:(max(i_group_tags(i_experiment_reps == er))*2)+.5);
        set(gca,'XTickLabels',geno_list{set_token}.colheaders); % X Labels
    end
    ylabel({'? Compressibility' ; '(per 500 modules)'},'Fontsize',32); % Y Labels
%     axis([0.5 (max(i_group_tags(i_experiment_reps == er))*2)+.5 ...
%         (min(scrap(1,:)) - (min(scrap(1,:))*0.05)) (max(scrap(2,:)) + (max(scrap(2,:))*0.05))]);
      axis([0.5 (max(i_group_tags(i_experiment_reps == er))*2)+.5 -0.022 0.0415]); % hard coded 
end

clear er set_token g scrap counter data

%% Compressibility Two Way ANOVA
dn_hour(1:14) = 1; dn_hour(15:24) = 2; dn_hour(25:38) = 3; dn_hour(39:48) = 4; 

for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % settings
    clear data; 
    
    for t = 1:max(dn_hour) % for each day/night 
        data(:,t) = nanmean(compressibility(i_experiment_reps == er,dn_hour == t),2); 
    end % average for each fish 
    data(:,find(sum(isnan(data)) == size(data,1),1,'first'):end) = []; % remove nans on the end 

    % Grouping Variables
    anova_group = repmat(i_group_tags(i_experiment_reps==er),...
        [size(data,2),1])'; % groups
    
    anova_experiment = repmat(i_experiment_tags(i_experiment_reps==er),...
        [size(data,2),1])'; % experiments
    
    anova_time = [];
    for t = time_window{set_token}(1):time_window{set_token}(2) % For each time window
        anova_time = [anova_time ; ones(sum(i_experiment_reps==er),1)*mod(t,2)];
        % Allocate alternating zeros and ones to each time window
    end
    anova_time = anova_time';
    
    % Development Grouping Variable
    if size(days_crop{set_token}(days{set_token}),2) == ...
            size(nights_crop{set_token}(nights{set_token}),2) ...
            && size(days_crop{set_token}(days{set_token}),2) > 1 % If there are an equal number of windows (>1)
        
        anova_development = []; % development
        anova_development = zeros(1,size(anova_group,2)); % Pre-allocate
        d = 1:size(anova_development,2)/(size(time_window{set_token}(1):...
            time_window{set_token}(2),2)/2):...
            size(anova_development,2); % divide into "24h" windows
        for t = 1:size(d,2)-1
            anova_development(d(t):d(t+1)-1) = t;
        end
    else
        anova_development = ones(size(anova_experiment)); % use all ones
    end
    
    % Comparison
    data = data(:)';
    
    [twa.rc.p{er}(:,1),~,twa.rc.stats{er}] = anovan(data,...
        {anova_group,anova_time,anova_development,anova_experiment},...
        'display','off','model','full');
    
end

clear er anova_group anova_experiment anova_time t anova_development d data

%% Module Classifiers 
    % Adapted from mRMR_Ideas.m
    
% Structure Data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'bout_proportions'); % fish x modules x days/nights
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'i_experiment_reps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'i_group_tags');

bout_proportions{1,1} = permute(bout_proportions{1,1},[1 3 2]); % now: fish x days/nights x modules
bout_proportions{2,1} = permute(bout_proportions{2,1},[1 3 2]);

%% WT Day and Night 
er = 1; 
mRMR_data{er,1} = double([reshape(bout_proportions{2,1}(i_experiment_reps == er,3:6,:),[],5,1) ...
    reshape(bout_proportions{1,1}(i_experiment_reps == er,3:6,:),[],5,1)]);
mRMR_tw{er,1} = [ones(sum(i_experiment_reps == 1),1) ; ones(sum(i_experiment_reps == 1),1)*2]; 
mRMR_tw{er,1} = [mRMR_tw{er,1} ; mRMR_tw{er,1}]; 

%% WT Day 5 vs Day 6
er = 1; 
mRMR_data{er,1} = double([reshape(bout_proportions{2,1}(i_experiment_reps == er,[3 5],:),[],5,1) ...
    reshape(bout_proportions{1,1}(i_experiment_reps == er,[3 5],:),[],5,1)]);
mRMR_tw{er,1} = [ones(sum(i_experiment_reps == 1),1) ; ones(sum(i_experiment_reps == 1),1)*2]; 

%% WT Night 5 vs Night 6 
er = 1; 
mRMR_data{er,1} = double([reshape(bout_proportions{2,1}(i_experiment_reps == er,[4 6],:),[],5,1) ...
    reshape(bout_proportions{1,1}(i_experiment_reps == er,[4 6],:),[],5,1)]);
mRMR_tw{er,1} = [ones(sum(i_experiment_reps == 1),1) ; ones(sum(i_experiment_reps == 1),1)*2]; 

%% Hcrtr Day and Night (Pairwise)
er = 2; 
mRMR_data{er,1} = double([reshape(bout_proportions{2,1}(i_experiment_reps == er,3:6,:),[],5,1) ...
    reshape(bout_proportions{1,1}(i_experiment_reps == er,3:6,:),[],5,1)]);
mRMR_tw{er,1} = repmat(i_group_tags(i_experiment_reps == er),4,1);
comps = 10; 

%% Hcrtr Day (Pairwise)
er = 2; 
mRMR_data{er,1} = double([reshape(bout_proportions{2,1}(i_experiment_reps == er,[3 5],:),[],5,1) ...
    reshape(bout_proportions{1,1}(i_experiment_reps == er,[3 5],:),[],5,1)]);
mRMR_tw{er,1} = repmat(i_group_tags(i_experiment_reps == er),2,1);
comps = 10; 

%% Hcrtr Night (Pairwise) 
er = 2; 
mRMR_data{er,1} = double([reshape(bout_proportions{2,1}(i_experiment_reps == er,[4 6],:),[],5,1) ...
    reshape(bout_proportions{1,1}(i_experiment_reps == er,[4 6],:),[],5,1)]);
mRMR_tw{er,1} = repmat(i_group_tags(i_experiment_reps == er),2,1);
comps = 10; 

%% Melatonin (Day) 
er = 3; 
mRMR_data{er,1} = double([reshape(bout_proportions{2,1}(i_experiment_reps == er,1,:),[],5,1) ...
    reshape(bout_proportions{1,1}(i_experiment_reps == er,1,:),[],5,1)]);
mRMR_tw{er,1} = i_group_tags(i_experiment_reps == er); 

%% PTZ (Day) 
er = 4; 
mRMR_data{er,1} = double([reshape(bout_proportions{2,1}(i_experiment_reps == er,1,:),[],5,1) ...
    reshape(bout_proportions{1,1}(i_experiment_reps == er,1,:),[],5,1)]);
mRMR_tw{er,1} = i_group_tags(i_experiment_reps == er); 

%% "Leave one out" Classifiers 

comps = 10; % number of modules to use (all)

tic
counter = 1;
for g_one = min(mRMR_tw{er,1}):max(mRMR_tw{er,1}) % for each group
    
    tags = ones(size(mRMR_tw{er,1}))*2; % all data (marked as 2)
    tags(mRMR_tw{er,1} == g_one) = 1; % group of interest (marked as 1) 
    
    % mRMR
    [comps_v{er,1}(counter,:)] = mrmr_miq_d(...
        zscore(mRMR_data{er,1}),tags,comps);
    
    % Classifiers
    for s = 1:comps % for each comp sequence
        % Fit a linear classifier as you add features
        % Using 10 fold cross validation
        % Hold 10% of mRMR_data back by default
        Mdl = fitcdiscr(...
            zscore(mRMR_data{er,1}(:,...
            comps_v{er,1}(counter,1:s))),...
            tags,...
            'DiscrimType','linear','CrossVal','on');
        Mdl_loss{er,1}(counter,s) = kfoldLoss(Mdl);
        Mdl_loss{er,2}(counter,s) = nanstd(kfoldLoss(Mdl,'Mode','individual'));
    end
    
    %Minimal Feature Space
    mRMR_ms(er,counter) = find(smooth(Mdl_loss{er,1}(counter,:),3) == ...
        min(smooth(Mdl_loss{er,1}(counter,:),3)),1,'first');
    
    if counter == 1 % remove WT data 
       mRMR_data{er,1}(mRMR_tw{er,1} == 1,:) = [];
       mRMR_tw{er,1}(mRMR_tw{er,1} == 1,:) = []; 
    end
    
    disp(num2str(counter));
    counter = counter + 1;
end

toc

%% Module Classifiers - Load Data (Experiment Colours) 
load('C:\Users\Marcus\Documents\Thesis\Tables\Draft_1_Mdl_Loss\Mdl_Loss_Data.mat'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\mRMR_Comparisons\Modules\module_data.mat');
motif_data = data([1:3 30:end],1:2);
tags = MdlLossTags([1:3 30:end],1);
clear data MdlLossTags; 

% Colormaps 
cmap(1:3,:) = repmat(([1 1 1]*(1-(1/(9)^.5))),3,1); % grey
cmap(4:12,:) = repmat([1 0.5 0],9,1); % orange   
cmap(13:19,:) = repmat([0 0.6178 0.8519],7,1); % blue 
cmap(20:23,:) = repmat([0.6863 0.2078 0.2784],4,1); % red

%% Module Classifiers - Figure
  
figure;
hold on; set(gca,'FontName','Calibri');
plot([0 54],[0 54],'--k','linewidth',3); % Hard coded line 

for i = 1:length(tags)
    errorbar(module_data(i,1), motif_data(i,1),...
        motif_data(i,2)/2,motif_data(i,2)/2,...
        module_data(i,2)/2,module_data(i,2)/2,...
        'color',cmap(i,:),...
        'marker','o','linewidth',3);
end

box off; set(gca,'Layer','top'); set(gca,'Fontsize',32); % Format
xlabel({'Module' ; 'Classification Error (%)'},'Fontsize',32); % Y Labels
ylabel({'Motif' ; 'Classification Error (%)'},'Fontsize',32); % Y Labels
axis([0 54 0 54]); 


%% Module Classifiers - Load Data (Group Colours) 
load('C:\Users\Marcus\Documents\Thesis\Tables\Draft_1_Mdl_Loss\Mdl_Loss_Data.mat'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\mRMR_Comparisons\Modules\module_data.mat');
motif_data = data([1:3 30:end],1:2);
tags = MdlLossTags([1:3 30:end],1);
clear data MdlLossTags; 

% Colormaps 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'cmap'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'cmap_2');
cmap{4} = [cmap{4}(2,:) ; cmap{4}(3,:) ; [1 0.5 0]]; % het, hom, orange 
cmap = [cmap{1} ; cmap_2{1} ; cmap{4} ; cmap{4} ; cmap{4} ; ...
    cmap{6} ; cmap{7}];
experiment_reps = zeros(length(cmap),1); 
experiment_reps(1:3) = 1; 
experiment_reps(4:6) = 2; 
experiment_reps(13:19) = 3; 
experiment_reps(20:end) = 4;

%% Module Classifiers - Figure
  
figure;
for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % settings
    subplot(1,4,er); counter = 1; % counts groups for plots
    hold on; set(gca,'FontName','Calibri'); clear scrap;
    
    for i = find(experiment_reps == er)'
        errorbar(module_data(i,1), motif_data(i,1),...
            motif_data(i,2)/2,motif_data(i,2)/2,...
            module_data(i,2)/2,module_data(i,2)/2,...
            'color',cmap(i,:),...
            'marker','o','linewidth',3);
    end
        
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
    ylabel({'Module' ; 'Classification Error (%)'},'Fontsize',32); % Y Labels
    ylabel({'Motif' ; 'Classification Error (%)'},'Fontsize',32); % Y Labels
    axis tight
    %axis([0.5 (max(i_group_tags(i_experiment_reps == er))*2)+.5 -0.022 0.0415]); % hard coded
end

clear er set_token g scrap counter data

%% PCA Space  
load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\180515_1.mat');
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'idx_numComp_sorted');
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'cmap_cluster');

options = statset('MaxIter',max_its); % Max number of GMM iterations

for c = 1:max(idx_numComp_sorted{1,1}) % for each module
    
    GMModels{c} = fitgmdist(X{1,1}(idx_numComp_sorted{1,1} == c,1:2),1,...
        'Options',options,'RegularizationValue',...
        rv,'Replicates',GMM_reps); % Fit k gaussians to sample, GMM_reps times
    
    GMModles_pdf{c} = pdf(GMModels{c},X{1,1}(idx_numComp_sorted{1,1} == c,1:2));
end

%% Figure

for c = max(idx_numComp_sorted{1,1}):-1:1 % for each module
    scatter3(X{1,1}(idx_numComp_sorted{1,1} == c,1),...
        X{1,1}(idx_numComp_sorted{1,1} == c,2),...
        GMModles_pdf{c},2,'filled',...
        'MarkerEdgeColor',cmap_cluster{1,1}(c,:),...
        'MarkerFaceColor',cmap_cluster{1,1}(c,:)); 
    
    hold on;
end

set(gca,'FontName','Calibri');
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
xlabel('PC1','Fontsize',32); % x Label
ylabel('PC2','Fontsize',32); % Y Label
zlabel('Probability','Fontsize',32);
axis tight

%% Bout Shapes Figure - Edit 

% Load Data
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180522.mat',...
    'bouts','cmap_cluster','fps','numComp','unit_conversion',...
    'parameter_dists','units');
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'cells');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'idx_numComp_sorted');

% Cropping Pdfs by Module Minimums and Maximums 
scrap(:,1) = grpstats(cells{2,1}(:,3),idx_numComp_sorted{2,1},'min'); 
scrap(:,2) = grpstats(cells{2,1}(:,3),idx_numComp_sorted{2,1},'max'); 
scrap(1,:) = []; 

for c = 1:size(scrap,1) % for each module 
    parameter_dists{2,1}(c,1:scrap(c,1)-2) = NaN;
    parameter_dists{2,1}(c,scrap(c,2)+2:end) = NaN;
end 

% Active Shapes 
% Plotting variables
for b = 1:size(bouts,2) % for each cluster
    bs_l(b) = size(bouts{1,b},2); % find the max length
end

bs_l = max(bs_l); % max length

figure;
subplot(1,2,1); hold on; set(gca,'FontName','Calibri');
for b = 1:size(bouts,2) % for each bout type 
    legend_lines(b) = shadedErrorBar(1:(size(bouts{1,b},2)+2),...
        nanmean([zeros(size(bouts{1,b},1),1) bouts{1,b} zeros(size(bouts{1,b},1),1)]),...
        nanstd([zeros(size(bouts{1,b},1),1) bouts{1,b} zeros(size(bouts{1,b},1),1)])/sqrt(size(bouts{1,1},1)),...
        'lineProps',{'Color',cmap_cluster{1,1}(b,:),'LineWidth',3});
    legend_cols(b) = legend_lines(b).mainLine; % Store color
end

box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
axis([1 bs_l ylim]);
set(gca,'XTick',2:2:bs_l);
set(gca,'XTickLabels',{(round((1:2:bs_l)/fps{1},2,'decimals'))}); % hard coded
xlabel('Time (Seconds)','Fontsize',32);
ylabel('? Pixels','Fontsize',32);
title(horzcat('Active Modules'),'Fontsize',32);

% Legend 
[~,icons,plots,~] = legend(flip(legend_cols),string(numComp(1):-1:1),'Location','best');
legend('boxoff'); 
set(icons(1:numComp(1)),'Fontsize',32) ; set(plots,'LineWidth',3);

% Inactive Module Fits 
clear legend_lines legend_cols; 
subplot(1,2,2); hold on;
box off; set(gca,'Layer','top'); set(gca,'Fontsize',32); set(gca,'FontName','Calibri'); % Set Font

for c = size(scrap,1):-1:1 % for each module 
    legend_lines(c) = plot(parameter_dists{2,1}(c,:),...
        'Color',cmap_cluster{2,1}(c,:),'LineWidth',0.001); % draw a thin line for the legned 
    area(parameter_dists{2,1}(c,:),...
         'EdgeColor',[1 1 1],'FaceColor',cmap_cluster{2,1}(c,:)); 
end 

set(gca,'XTick',[1 10 fps{1}*60]); 
set(gca,'XTickLabel',{0.04 0.40 60.00})
axis([1 fps{1}*60 0 max(parameter_dists{2,1}(:))]);
set(gca,'XScale','log'); % set log axis
xlabel('Time (Seconds)','Fontsize',32);
ylabel('Probability','Fontsize',32); % Y label

% Legend 
[~,icons,~,~] = legend(flip(legend_lines),string(numComp(2):-1:1),'Location','best');
legend('boxoff'); 
set(icons(1:numComp(2)),'Fontsize',32) ; 
set(icons(6:end),'LineWidth',3);
title(horzcat('Inactive Modules'),'Fontsize',32);

clear b bs_l legend_lines legend_cols icons plots s p crop k  

%% Video 

v = VideoReader('D:\Behaviour\SleepWake\Videos\171013_18_19\171013_18_19_c97_0001.avi');

vidWidth = v.Width;
vidHeight = v.Height;

mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint64'),...
    'colormap',[]); % pre-allocate

k = 1;
while hasFrame(v)
    mov(k).cdata = readFrame(v);
    k = k+1;
end

%% Delta Px Data 

load('D:\Behaviour\SleepWake\Videos\171013_18_19\delta_px_sq.mat');

% Normalise Data
delta_px_sq = delta_px_sq./max(delta_px_sq); % normalise 
delta_px_sq(isnan(delta_px_sq) == 1) = 0; % remove NaN Values 

% Rescale 
delta_px_sq = delta_px_sq * 50; % hard coded scalar 

% Remove Extra Frames (hard coded)
delta_px_sq(1:5,:) = [];

fps = v.framerate; % get framerate

%% Video
%Positions Vector
sc = [90 , 150]; % Starting coordinates (hard coded)
x = 72; % X distance between wells (hard coded)
y = 74; % Y distance between wells (hard coded)

p = 1; % For each position
for r = 1:8 % For each row (Hard coded)
    for c = 1:12 % For each column (Hard coded)
        well_positions(p,:) = [sc(1)+((c*x)-x), sc(2)+((r*y)-y)];
        p = p+1;
    end
end

% Visually Check Grid 
% figure; 
% imagesc(mov(1).cdata);
% hold on
% scatter(well_positions(:,1),well_positions(:,2));

% Making the Video
data_frames = fps:(fps*30 + fps - 1); % Take 30s of video

% Pre-allocate a strucutre to store video frames
s(length(data_frames)) = struct('cdata',[],'colormap',[]);

set(0,'DefaultFigureVisible','off'); % Suppress figures from popping up

for i = data_frames %For each frame
    
    h = axes('position',[0 0 1 1]); % Initial axis
    image(mov(i).cdata); hold on; %Draw the frame
    set(gca,'XTick',[]); % Remove x-axis ticks
    set(gca,'YTick',[]); % Remove y-axis ticks
    
    % Fish Behaviour
    for f = 1:size(delta_px_sq,2) % for each fish
        plot(well_positions(f,1):(well_positions(f,1)+fps-1),...
            (well_positions(f,2) - delta_px_sq((i-fps+1):i,f)),...
            'color',[0.3922 0.5843 0.9294],'linewidth',3); % Plot 1s of data before this frame
        
        % To this frame
        scatter(well_positions(f,1)+fps-1,well_positions(f,2) - delta_px_sq(i,f),...
            20,'MarkerFaceColor',...
            'k','MarkerEdgeColor','none'); % Draw a black Dot at the current value
    end
    
    % Storing frames with drawings on top
    drawnow;
    s(i+1-fps) = getframe(h); % Grab the frame 
    disp(num2str(i+1-fps)); % Report back the frames
    close all % ! % Note this two commands incrase the writing speed
    hold off % !  % By decreasing the amount of data held in RAM
    
end

set(0,'DefaultFigureVisible','on') % Allow figures to pop up again

figure;
% Open a new figure and play the movie from the structure
movie(s,1,v.FrameRate) % Play the movie at the native framerate

%Write the Video
vOut = VideoWriter('C:\Users\Marcus\Desktop\Video_Test.avi','MPEG-4');
vOut.FrameRate = v.FrameRate;
vOut.Quality = 100; 
%vOut.LosslessCompression = true; 
open(vOut)
for k = 1:numel(s)
    writeVideo(vOut,s(k));
end
close(vOut)

%% Cluster Parameters - Fits Figure (Inactive Bout Length Corrected) 

% Load data 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat',...
    'cells','parameters','unit_conversion','parameter_dists',...
    'cmap_cluster','units','numComp');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'idx_numComp_sorted');

% Cropping Pdfs by Module Minimums and Maximums 
scrap(:,1) = grpstats(cells{2,1}(:,3),idx_numComp_sorted{2,1},'min'); 
scrap(:,2) = grpstats(cells{2,1}(:,3),idx_numComp_sorted{2,1},'max'); 
scrap(1,:) = []; 

for c = 1:size(scrap,1) % for each module 
    parameter_dists{2,1}(c,1:scrap(c,1)-2) = NaN;
    parameter_dists{2,1}(c,scrap(c,2)+2:end) = NaN;
end 
 
%% Figure 
figure;
counter = 1; % start a counter 
for s = 1:2 % for active & inactive
    for p = 3:size(cells{s,1},2) % For each parameter
        % Figure Settings
        subplot(2,4,counter); hold on;
        box off; set(gca,'Layer','top'); set(gca,'Fontsize',12); set(gca,'FontName','Calibri'); % Set Font
        if s ~= 2 % for active parameters 
            title(parameters{1}{p-2}); % Add title
        else % for inactive bout length 
            title(parameters{1}{10}); % Add title   
        end
        
        % Plot 
        if min(cells{s,1}(:,p)) == 0 % shift from zero (for log scaling)
            cells{s,1}(:,p) = cells{s,1}(:,p) + 1;
            off_set = 1;
        end
        
        crop = max(cells{s,1}(:,p)); 
        
        if s == 1 % if active 
            for k = 1:numComp(s) % for each cluster
                plot((min(cells{s,1}(:,p)):crop)/unit_conversion{1}(s,p-2),...
                    parameter_dists{s,p-2}(k,:),'color',cmap_cluster{s,1}(k,:),'linewidth',3);
            end
        else % if inactive 
            for k = numComp(s):-1:1 % for each cluster
                plot((min(cells{s,1}(:,p)):crop)/unit_conversion{1}(s,p-2),...
                    parameter_dists{s,p-2}(k,:),'color',cmap_cluster{s,1}(k,:),'linewidth',0.001);
                area((min(cells{s,1}(:,p)):crop)/unit_conversion{1}(s,p-2),...
                    parameter_dists{s,p-2}(k,:),...
                    'EdgeColor',[1 1 1],'FaceColor',cmap_cluster{2,1}(k,:));
            end
        end 

        % Axes 
        try 
            set(gca,'XTick',...
                [min(cells{s,1}(:,p))/unit_conversion{1}(s,p-2), (min(cells{s,1}(:,p))/unit_conversion{1}(s,p-2))*10,...
                crop/unit_conversion{1}(s,p-2)]); % set x tick labels
        catch 
            set(gca,'XTick',...
                [min(cells{s,1}(:,p))/unit_conversion{1}(s,p-2), (crop/unit_conversion{1}(s,p-2))/2,...
                crop/unit_conversion{1}(s,p-2)]); % set x tick labels
        end 
        axis([min(cells{s,1}(:,p))/unit_conversion{1}(s,p-2) crop/unit_conversion{1}(s,p-2) ...
            min(parameter_dists{s,p-2}(:)) max(parameter_dists{s,p-2}(:))]); % Set axis limits
        
        % Set decimal places depending on units
        if unit_conversion{1}(s,p-2) > 1
            xtickformat('%.2f');
        else
            xtickformat('%.0f');
        end
        
        if exist('off_set','var') == 1  % check to see if data was shited from zero
            cells{s,1}(:,p) = cells{s,1}(:,p) - 1; % remove 1 from data
            set(gca,'XTickLabels',string(cellfun(@str2num,get(gca,'XTickLabels')) - 1)); % remove 1 from x-labels
            clear off_set;
        end
  
        set(gca,'XScale','log'); % set log axis
        xlabel(units{1}(p-2),'Fontsize',12); % X labels
        ylabel('Probability','Fontsize',12); % Y label
        
        counter = counter + 1; % add to counter 
        
        % Legend Plot 
        if counter == 8 
            subplot(2,4,counter); hold on;
            box off; set(gca,'Layer','top'); set(gca,'Fontsize',12); set(gca,'FontName','Calibri'); % Set Font
            title('Legend'); 
            
            for s_2 = 1:2 % for active & inactive 
                for k = 1:numComp(s_2) % for each cluster 
                    scatter(k,mod(s_2,2),300,...
                        'markerfacecolor',cmap_cluster{s_2,1}(k,:),...
                        'markeredgecolor',cmap_cluster{s_2,1}(k,:));
                end 
            end 
            
            axis([0 max(numComp)+1 -1 2]);
            xlabel('Module Number','Fontsize',12); % X labels
            set(gca,'XTick',1:max(numComp)); 
            set(gca,'YTick',[0 1]); 
            set(gca,'YTickLabels',{'Inactive','Active'},'Fontsize',12);
        end 
        
    end
end

clear counter s p crop k s_2

%% Module Video - Preparing Model 

load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\Pre.mat', 'wake_cells');
load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\Pre.mat', 'fish_tags');

% Active 
tic
X{1,1} = []; % empty X
% Z-score each fishes data
for f = 1:max(fish_tags{1,1}) % for each fish
    X{1,1} = [X{1,1} ; zscore(wake_cells(fish_tags{1,1} == f,3:end))];
    if mod(f,100) == 0 % report every 100 fish 
        disp(horzcat('Completed ',num2str(f),' fish of ',...
            num2str(max(fish_tags{1,1}))));
    end
end
toc

mu = mean(X{1,1}); % mean Model feature values 

clearvars -except mu

%% Module Video - Load Clustered Data  
load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\Pre.mat','wake_cells','sleep_cells'); 
active = load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\180515_1.mat','idx');
inactive = load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\180515_2.mat','idx'); 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'numComp')
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'sleep_cells_nan_track')

% merge variables 
idx{1,1} = active.idx; idx{2,1} = inactive.idx; 

clear active inactive

% Remove NaN Values from Sleep matracies  
sleep_cells(sleep_cells_nan_track,3) = NaN; % sleep cells 
idx{2,1}(sleep_cells_nan_track,1) = NaN; % cluster assignments  

% Merge Variables for ease of Looping 
cells{1,1} = wake_cells; 
cells{2,1} = sleep_cells; % cells 

%% Module Video - Sorting Clusters by Mean Length 
             
for s = 1:2 % for active & inactive
    mean_cluster_length = nan(1,numComp(s),'single'); % pre-allocate
    for c = 1:numComp(s) % For each cluster
        mean_cluster_length(c) = nanmean(cells{s,1}(idx{s,1}==c,3));
        % Calculate mean bout length
    end
    
    [~,O] = sort(mean_cluster_length); % Sort by increasing bout length
    
    cluster_order(s,:) = O; 
    
    clear s mean_cluster_length c O;
end

clearvars -except mu cluster_order  

%% Module Video - Load Video  

v = VideoReader('D:\Behaviour\SleepWake\Videos\171013_18_19\171013_18_19_c97_0001.avi');

vidWidth = v.Width;
vidHeight = v.Height;

mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint64'),...
    'colormap',[]); % pre-allocate

k = 1;
while hasFrame(v)
    mov(k).cdata = readFrame(v);
    k = k+1;
end

%% Module - Delta Px Data 

load('D:\Behaviour\SleepWake\Videos\171013_18_19\delta_px_sq.mat');

% Normalise Data
delta_px_sq = delta_px_sq./max(delta_px_sq); % normalise 
delta_px_sq(isnan(delta_px_sq) == 1) = 0; % remove NaN Values 

% Rescale 
delta_px_sq = delta_px_sq * 50; % hard coded scalar 

% Remove Extra Frames (hard coded)
delta_px_sq(1:5,:) = [];

fps = v.framerate; % get framerate

%% Modules - Extract Bouts 

% Pre-allocate 
wake_cells = cell(1,size(delta_px_sq,2)); % Wake Cells (bout parameters)  
sleep_cells = cell(1,size(delta_px_sq,2)); % Sleep Cells (bout parameters) 

% Finding transitions
delta_px_sq_scrap = delta_px_sq; 
delta_px_sq_scrap(delta_px_sq_scrap > 0) = 1; % Find active frames 
delta_px_sq_scrap = diff(delta_px_sq_scrap); % Diff to find transitions  
    % 1 = inactive to active 
    % -1 = active to inactive 
 
for f = 1:size(delta_px_sq,2) % For each fish 
    
    % Starts - ensures no bouts are lost at the start  
    if  delta_px_sq(1,f) > 0 % If active in first bin  
        wake_cells{1,f}(:,1) = [1 ; find(delta_px_sq_scrap(:,f) == 1)+1]; % Find active bout starts
        sleep_cells{1,f}(:,1) = (find(delta_px_sq_scrap(:,f) == -1)+1); % Find sleep bout starts 
    else % Ie. if inactive in first bin 
        wake_cells{1,f}(:,1) = find(delta_px_sq_scrap(:,f) == 1)+1; % Find active bout starts 
        sleep_cells{1,f}(:,1) = [1 ; (find(delta_px_sq_scrap(:,f) == -1)+1)]; % Find sleep bout starts 
    end 
    
    % Ends - ensures no bouts are lost at the end 
    if delta_px_sq(size(delta_px_sq,1),f) > 0 % If active in last bin 
        wake_cells{1,f}(:,2) = [find(delta_px_sq_scrap(:,f) == - 1);...
            size(delta_px_sq,1)]; % Find active bout ends
        sleep_cells{1,f}(:,2) = find(delta_px_sq_scrap(:,f) == 1); % Find sleep bout ends
    else 
        wake_cells{1,f}(:,2) = find(delta_px_sq_scrap(:,f) == - 1); 
        sleep_cells{1,f}(:,2) = [(find(delta_px_sq_scrap(:,f) == 1)) ; size(delta_px_sq,1)]; % Find sleep bout ends
    end
    
    % Parameter extraction 
    wake_cells{1,f}(:,3:8) = NaN; % Pre-allocate 
    wake_cells{1,f}(:,3) = (wake_cells{1,f}(:,2)+1) - wake_cells{1,f}(:,1); % Wake Bout Length 
    sleep_cells{1,f}(:,3) = (sleep_cells{1,f}(:,2)+1) - sleep_cells{1,f}(:,1); % Sleep Bout Length 

    % Active bouts 
    for b = 1:size(wake_cells{1,f},1) % For each active bout
        wake_cells{1,f}(b,4) = nanmean(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Mean
        wake_cells{1,f}(b,5) = nanstd(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Std
        wake_cells{1,f}(b,6) = nansum(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Total
        wake_cells{1,f}(b,7) = nanmin(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Min
        wake_cells{1,f}(b,8) = nanmax(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Max
    end
end 

% Concatenate variables
% cells{1,1} = []; 
% cells{2,1} = []; 
fish_tags{1,1} = []; 
fish_tags{2,1} = []; 
counter = 1;
for i = 1:size(wake_cells,2) % For each fish
%     cells{1,1} = [cells{1,1} ; wake_cells{1,i}]; % wake cells
%     cells{2,1} = [cells{2,1} ; sleep_cells{1,i}]; % sleep cells
    fish_tags{1,1} = [fish_tags{1,1} ; ones(size(wake_cells{1,i},1),1)*counter]; % wake fish tags
    fish_tags{2,1} = [fish_tags{2,1} ; ones(size(sleep_cells{1,i},1),1)*counter]; % sleep fish tags
    counter = counter + 1; % add to fish counter
end

wake_cells = cell2mat(wake_cells'); % check arrangement
sleep_cells = cell2mat(sleep_cells'); % check arrangmenet 

%% Preparing New Data   

% Active 
X{1,1} = []; % empty X
% Z-score each fishes data
for f = 1:max(fish_tags{1,1}) % for each fish
    X{1,1} = [X{1,1} ; zscore(wake_cells(fish_tags{1,1} == f,3:end))];
end

% Load Model 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'coeff');
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'score');
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'knee_dim'); 

% Fit To Model's PCA Basis 
% https://stackoverflow.com/questions/13303300/how-to-project-a-new-point-to-pca-new-basis#
X{1,1} = bsxfun(@minus,X{1,1}, mu); % subtract old means from new data 
X{1,1} = X{1,1}*coeff; % project new data into pca basis 
X{1,1} = X{1,1}(:,1:knee_dim); % crop to number of PC's 

% Inactive 
X{2,1} = []; % empty X
X{2,1} = sleep_cells(:,3); 

% Convert from frames to seconds  
X{2,1} = X{2,1}/fps; 

%% Assigning New Data to Clusters 

% Load variables 
nn = 50; 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'sample_a_n');
Mdl = load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'X'); 
Mdl = Mdl.X; % old data in PCA / Frame space 
Mdl{2,1} = Mdl{2,1}/fps; % convert from Frames to seconds  

for s = 1:2 % for active and inactive
    disp('Fitting data into clusters'); % report progress
    
    ea_idx_n = sample_a_n{s,1}(:,2); % normalised cluster indicies
    
    idx{s,1} = zeros(size(X{s,1},1),1,'single'); % allocate
    cks = [1:round(size(X{s,1},1)/1000):size(X{s,1},1) (size(X{s,1},1)+1)]; % break into 1000 chunks
    for i = 1:(length(cks) - 1) % for each chunk
        kn = knnsearch(Mdl{s,1}(sample_a_n{s,1}(:,1),:),X{s,1}(cks(i):(cks(i+1)-1),:),...
            'K',nn); % find nn nearest neighbours
        idx{s,1}(cks(i):(cks(i+1)-1),1) = ...
            mode(reshape(ea_idx_n(kn(:)),size(kn)),2); % assign to mode neighbour cluster
    end
end

% Convert back to frames
X{2,1} = X{2,1}*25; 

Mdl{2,1} = Mdl{2,1}*25; % convert from Frames to seconds  

%% Module - Video Sorting Clusters by Model Length 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'numComp'); 

for s = 1:2 % for active & inactive
    
    O = cluster_order(s,:);
    
    idx_numComp_sorted{s,1} = nan(size(idx{s,1},1),1,'single'); % Pre-allocate
    
    for c = 1:numComp(s)  % For each cluster
        idx_numComp_sorted{s,1}(idx{s,1} == O(c),:) = c;
        % Re-assign cluster numbers
    end
    
    clear s O c
end

clear idx

%% Module Video - Threading 

% Pre-allocation
threads = cell(max(fish_tags{1,1}),2); % fish x (clusters,times (start & stop)
idx_numComp_sorted{1,1} = idx_numComp_sorted{1,1} + max(idx_numComp_sorted{2,1}); % Assign higher numbers to the wake bouts
idx_numComp_sorted{2,1}(isnan(idx_numComp_sorted{2,1})) = 0; % replace nan-values with zero

for f = 1:max(fish_tags{1,1}) % For each fish
    % Pre-allocate
    % Data
    threads{f,1,1} = nan(size(find(fish_tags{1,1} == f),1) + ...
        size(find(fish_tags{2,1} == f),1),1,'single'); % clusters
    threads{f,2,1} = nan(size(threads{f,1,1},1),2,'single'); % times (start & stop)
    
    % Deterime starting state (a = active, b = inactive)
    if wake_cells(find(fish_tags{1,1} == f,1,'first'),1) == 1 % If the fish starts active
        a = 1; b = 2; % start active
    else % If the fish starts inactive
        a = 2; b = 1; % start inactive
    end
    
    % Fill in data
    % Clusters
    threads{f,1,1}(a:2:end,1) = idx_numComp_sorted{1,1}...
        (fish_tags{1,1} == f,1); % Fill in active clusters
    threads{f,1,1}(b:2:end,1) = idx_numComp_sorted{2,1}...
        (fish_tags{2,1} == f,1); % Fill in inactive clusters
    % Times
    threads{f,2,1}(a:2:end,1:2) = wake_cells...
        (fish_tags{1,1} == f,1:2); % Fill in active times
    threads{f,2,1}(b:2:end,1:2) = sleep_cells...
        (fish_tags{2,1} == f,1:2); % Fill in inactive times

end 

%% Module Video - States 

% Pre-allocate data
states = nan(size(delta_px_sq));

% Fill data in
for f = 1:size(states,2) % for each fish
    for b = 1:size(threads{f,1,1},1) % for each bout
        states(threads{f,2,1}(b,1):...
            threads{f,2,1}(b,2),f) = ...
            threads{f,1,1}(b,1); % Fill in cluster number
    end
end
    
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'cmap_cluster_merge');

%% Module Video - Video
%Positions Vector
sc = [90 , 150]; % Starting coordinates (hard coded)
x = 72; % X distance between wells (hard coded)
y = 74; % Y distance between wells (hard coded)

p = 1; % For each position
for r = 1:8 % For each row (Hard coded)
    for c = 1:12 % For each column (Hard coded)
        well_positions(p,:) = [sc(1)+((c*x)-x), sc(2)+((r*y)-y)];
        p = p+1;
    end
end

% Visually Check Grid 
% figure; 
% imagesc(mov(1).cdata);
% hold on
% scatter(well_positions(:,1),well_positions(:,2));

% Making the Video
data_frames = fps:(fps*30 + fps - 1); % Take 30s of video

% Pre-allocate a strucutre to store video frames
s(length(data_frames)) = struct('cdata',[],'colormap',[]);

set(0,'DefaultFigureVisible','off'); % Suppress figures from popping up

for i = data_frames %For each frame
    
    h = axes('position',[0 0 1 1]); % Initial axis
    image(mov(i).cdata); hold on; %Draw the frame
    set(gca,'XTick',[]); % Remove x-axis ticks
    set(gca,'YTick',[]); % Remove y-axis ticks
    colormap(cmap_cluster_merge); % set colormap 
    
    % Fish Behaviour
    for f = 1:size(delta_px_sq,2) % for each fish
        % Older Option 
        seq(1,:) = well_positions(f,1):(well_positions(f,1)+fps-1); 
        seq(2,:) = (well_positions(f,2) - delta_px_sq((i-fps+1):i,f))'; 
        seq(3,:) = states((i-fps+1):i,f)'; 
        for t = 1:(length(seq)-1)
            % Plot 
            plot([seq(1,t),seq(1,t+1)],[seq(2,t),seq(2,t+1)],...)
                'color',cmap_cluster_merge(max([seq(3,t) seq(3,t+1)]),:),...)
                'linewidth',3)
        end
        scatter(well_positions(f,1)+fps-1,well_positions(f,2) - delta_px_sq(i,f),...
            20,0,'MarkerFaceColor',cmap_cluster_merge(states(i,f),:),...
            'MarkerEdgeColor','none'); % Draw a black Dot at the current value
     end
    
    % Storing frames with drawings on top
    drawnow;
    s(i+1-fps) = getframe(h); % Grab the frame 
    disp(num2str(i+1-fps)); % Report back the frames
    close all % ! % Note this two commands incrase the writing speed
    hold off % !  % By decreasing the amount of data held in RAM
    
end

set(0,'DefaultFigureVisible','on') % Allow figures to pop up again

figure;
% Open a new figure and play the movie from the structure
movie(s,1,v.FrameRate) % Play the movie at the native framerate

%Write the Video
vOut = VideoWriter('C:\Users\Marcus\Desktop\Video_Test.avi','MPEG-4');
vOut.FrameRate = v.FrameRate;
vOut.Quality = 100; 
%vOut.LosslessCompression = true; 
open(vOut)
for k = 1:numel(s)
    writeVideo(vOut,s(k));
end
close(vOut)

%% Scrap Code 
%             % Scatter
%             scatter([seq(1,t) seq(1,t)],[seq(2,t) seq(2,t)],...
%                 'markerfacecolor',cmap_cluster_merge(seq(3,t),:),...
%                 'markeredgecolor',cmap_cluster_merge(seq(3,t),:));
%         % Old Code 
%         for t = 1:length(seq) % for each module in the sequence
%             if seq(t) <= numComp(1) % for the inactive modules
%                 plot([a (a+ibl(seq(t)))],[b b],...
%                     'color',cmap_cluster_merge(seq(t),:),'linewidth',5); % plot
%                 a = a + ibl(seq(t)); % add to time
%             else % for the active modules
%                 plot(a:(a+length(nanmean(bouts{1,seq(t)-numComp(1)}))+1),...
%                     [b ((nanmean(bouts{1,seq(t)-numComp(1)})/28)+b) b],...
%                     'color',cmap_cluster_merge(seq(t),:),'linewidth',5); % plot
%                 a = a + length(nanmean(bouts{1,seq(t)-numComp(1)})) + 1; % add to time
%             end
%         end
%     
%         % One Option 
%         clinep((well_positions(f,1):(well_positions(f,1)+fps-1)),...
%             (well_positions(f,2) - delta_px_sq((i-fps+1):i,f))',...
%             zeros(1,fps),states((i-fps+1):i,f)',4);
        
        % To this frame
%         scatter(well_positions(f,1)+fps-1,well_positions(f,2) - delta_px_sq(i,f),...
%             60,0,'MarkerFaceColor',cmap_cluster_merge(states(i,f),:),...
%             'MarkerEdgeColor','none'); % Draw a black Dot at the current value
