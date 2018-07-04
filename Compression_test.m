%% Fine-er grain compression 
    % For the lighting transitions 
    
% Load data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'threads');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'idx_numComp_sorted');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'cmap');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'cmap_2');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'i_group_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'i_experiment_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'time_window');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'lb'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'i_experiment_reps');

% Remove starting transition from short experiments 
time_window{6} = 2; 
time_window{7} = 2; 

%% Chunk Compression 
% Compress cluster chunks 

% settings
step = 500; % chunk size 
sMax = max(idx_numComp_sorted{1,1}) + 1; % Maximum states + 1 (first new symbol)
nMax = 10; % Maximum n-grams to consider

% allocate
chunks = cell(size(threads,1),1); % fish x 1
totSavings_cells = cell(size(threads,1),1); % fish x 1

% memory saving 
threads = threads(:,:,1); % removed shuffled data (eases worker memory) 

tic
errors = []; 
parfor f = 1:size(threads,1) % for each fish
    counter = 1; % counts lb boundaries
    set_token = i_experiment_tags(f); % used for each experiments sets settings
    
    for bound = min(time_window{set_token}):max(time_window{set_token}) % for each transition
        l_b = find(threads{f,2,1}(:,1) < lb{set_token}(bound),1,'last');
        
        try
            chunks{f,1}(counter,:) = threads{f,1,1}(l_b-749:l_b+1250,1);
            
            for t = 1:(length(chunks{f,1})-step) % for each chunk
                [~,~, totSavings_cells{f,1}(counter,t)] = ...
                    compressSequenceNFast(chunks{f,1}(counter,t:(t+(step-1))),...
                    sMax,nMax); % compress returning totSavings
            end
            
        catch % catch animals with few modules before the transition 
            chunks{f,1}(counter,1:2000) = NaN;
            totSavings_cells{f,1}(counter,1:1500) = NaN; 
            errors = [errors ; horzcat(num2str(f),'f,',num2str(counter))]; 
        end
        
        counter = counter + 1;
        
    end
    
    disp(horzcat('Grain-Compressed Fish ',num2str(f)));
    
end
toc

%% Calculate Compressibility

for f = find(i_experiment_reps == er)'
    
    for l_b = 1:size(totSavings_cells{f,1},1) % for each bound
        compressibility(f,:,l_b) = totSavings_cells{f,1}(l_b,:)/step; 
    end 
    
end

%% Compressibility Figure 

er = 1;
%cmap_temp = [0 0 0 ; 0 0 0 ; 1 0.5 0 ; cmap_2{1}(2,:)];
cmap_temp = [0 0 0 ; 0 0 0 ; cmap_2{1}(1,:) ; cmap_2{1}(2,:)];

cmap_temp(1,:) = cmap_temp(3,:)+(1-cmap_temp(3,:))*(1-(1/(3)^.5)); 
cmap_temp(2,:) = cmap_temp(4,:)+(1-cmap_temp(4,:))*(1-(1/(3)^.5)); 

figure; hold on; 
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32); 

for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
    
    for l_b = 1:size(compressibility,3) % for each l_b boundary
        legend_lines(l_b) = shadedErrorBar(1:size(compressibility,2),...
            nanmean(compressibility(i_group_tags(i_experiment_reps == er) == g,:,l_b)),...
            nanstd(compressibility(i_group_tags(i_experiment_reps == er) == g,:,l_b))...
            /sqrt(sum(i_group_tags(i_experiment_reps == er) == g)),...
            'lineprops',{'color',cmap_temp(l_b,:),'linewidth',3});
        
        legend_cols(l_b) = legend_lines(l_b).mainLine; % Store color
        
    end
    
end
axis tight
l = plot([size(compressibility,2)/2 size(compressibility,2)/2],ylim,...
    'color',[1 0.5 0],'lineStyle','--','linewidth',3); 
uistack(l,'bottom'); 
l = plot([(size(compressibility,2)/2)-step (size(compressibility,2)/2)-step],ylim,...
    'color','k','lineStyle','--','linewidth',1.5); 
uistack(l,'bottom'); 
l = plot([(size(compressibility,2)/2)+step (size(compressibility,2)/2)+step],ylim,...
    'color','k','lineStyle','--','linewidth',1.5); 
uistack(l,'bottom'); 


ylabel({'Compressibility' ; '(per 500 modules)'},'Fontsize',32); % Y Labels
set(gca,'XTick',[(size(compressibility,2)/2)-step size(compressibility,2)/2 ... 
    (size(compressibility,2)/2)+step]); 
xlabel('Module','Fontsize',32); % ylabel
legend([legend_cols(1) legend_cols(3) legend_cols(2) legend_cols(4)],...
    {'D:L (5dpf)','D:L (6dpf)','L:D (5dpf)','L:D (6dpf)'},...
    'Fontsize',32,'Location','NorthEast')
legend('boxoff'); 
