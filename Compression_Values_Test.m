%% Compression Values Test 

% settings
step = 500;
sMax = 11; % Maximum states + 1 (first new symbol)
nMax = 10; % Maximum n-grams to consider

%% Shuffled 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'threads'); 

%% Chunk Compression 
% Compress cluster chunks 
    % Note that threads{f,3,1} now contains hours and NaN values 
    % outside of the time windows of interest   

% allocate
chunks = cell(size(threads,1),1); % fish x 1
totSavings_cells = cell(size(threads,1),1); % fish x 1

for f = 1:size(threads,1) % for each fish
        chunks{f,1} = find(threads{f,3,1} == 1,1,'first'):...
        step:find(isnan(threads{f,3,1}) == 0,1,'last'); % chunk data  
end 

tic
parfor f = 1:size(threads,1) % for each fish
    
    for tc = 1:size(threads,3) % for real and shuffled data
        for t = 1:(length(chunks{f,1})-1) % for each chunk
            [~,~, totSavings_cells{f,1}(tc,t)] = ...
                compressSequenceNFast(threads{f,1,tc}(chunks{f,1}(1,t):(chunks{f,1}(t+1)-1))',...
                sMax,nMax); % compress returning totSavings
        end
    end
    
    disp(num2str(f));
    
end
toc

clear f tc t

%% Load data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Compression_Values_Test.mat');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'experiment_reps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'i_group_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'i_experiment_reps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'cmap');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'geno_list');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180524_Hours.mat', 'i_experiment_tags');

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

%% Random Sampling (compressibility_r) 

er = 1; 
for f = 1:length(i_group_tags(i_experiment_reps == er)) % for an equal number of samples 
    for t = 1:2 % for fictive day/night data 
        a = datasample(1:2,1); % choose a starting state (active or inactive)
        b = setdiff(1:2,a); % " " 
        scrap = nan(1,step,'single'); % allocate 
        scrap(a:2:end) = datasample(1:5,step/2,'Replace',true); % fill inactive data
        scrap(b:2:end) = datasample(6:10,step/2,'Replace',true); % fill active data
        [~,~,compressibility_r(f,t)] = compressSequenceNFast(scrap,sMax,nMax); % compress 
    end
end

compressibility_r = compressibility_r/step; % calculate compression  

%% Compressibility Day vs Night
    % Note 180612: make sure that dn_hour is set to 1-2
    
figure;
for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % settings
    ax(er) = subplot(2,2,er); counter = 1; % counts groups for plots
    hold on; set(gca,'FontName','Calibri'); clear scrap;
    
    for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
        clear data;
        % average day per fish
        data(:,1) = nanmean(compressibility(i_experiment_reps == er & i_group_tags == g,dn_hour == 1,1),2);
        % average night per night
        data(:,2) = nanmean(compressibility(i_experiment_reps == er & i_group_tags == g,dn_hour == 2,1),2);
        
        if er == 1
            shuffled = [nanmean(nanmean(compressibility(i_experiment_reps == er & ...
            i_group_tags == g,dn_hour == 1,2:end),2),3) ... 
            nanmean(nanmean(compressibility(i_experiment_reps == er & ...
            i_group_tags == g,dn_hour == 2,2:end),2),3)]; 
        
            data = [compressibility_r shuffled data]; 
            
            % random data (compressibility_r)
            plot([counter,counter+1],data(:,1:2),...
                'color',[1 0.5 0] + (1-[1 0.5 0])*(1-(1/(5)^.5)),'linewidth',1.5);
             errorbar([counter,counter+1],nanmean(data(:,1:2)),nanstd(data(:,1:2)),...
                'color',[1 0.5 0],'linewidth',3);
            counter = counter + 2; 
            
            % shuffled data 
            plot([counter,counter+1],data(:,3:4),...
                'color',([1 1 1]*(1-(1/(5)^.5))),'linewidth',1.5);
            errorbar([counter,counter+1],nanmean(data(:,3:4)),nanstd(data(:,3:4)),...
                'color',[0 0 0],'linewidth',3); 
            counter = counter + 2;
            
            % real data 
            plot([counter,counter+1],data(:,5:6),...
                'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
            errorbar([counter,counter+1],nanmean(data(:,5:6)),nanstd(data(:,5:6)),...
                'color',cmap{set_token}(g,:),'linewidth',3);           
           
        else 
            plot([counter,counter+1],data,...
                'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
            errorbar([counter,counter+1],nanmean(data),nanstd(data),...
                'color',cmap{set_token}(g,:),'linewidth',3);
            counter = counter + 2; % add to counter
        end 

    end
    
    y_lims(1,er) = min(data(:));
    y_lims(2,er) = max(data(:));
    
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
    if er == 1 % for the WT Data
        set(gca, 'XTick', [1.5 3.5 5.5]); % set X-ticks
        set(gca,'XTickLabels',{'Random','Shuffled','Real'}); % X Labels
        set(ax(er),'XLim',[.5 counter+1.5]); % set x limits for all subplots
    else % for the other experiments
        set(gca,'XTick',1.5:2:(max(i_group_tags(i_experiment_reps == er))*2)+.5);
        set(gca,'XTickLabels',geno_list{set_token}.colheaders); % X Labels
        set(ax(er),'XLim',[.5 (max(i_group_tags(i_experiment_reps == er))*2)+.5]); % set x limits for all subplots
    end
    xlabel('Data','Fontsize',32);
    ylabel({'Compressibility' ; '(per 500 modules)'},'Fontsize',32); % Y Labels

end

set(ax,'Ylim',[min(y_lims(:)) max(y_lims(:))]); % set y limits for all subplots   

clear er set_token g scrap counter data

%% WT N-WAY ANOVA
dn_hour(1:14) = 1; dn_hour(15:24) = 2; dn_hour(25:38) = 3; dn_hour(39:48) = 4;

er = 1; g = 1; 
set_token = find(experiment_reps == er,1,'first'); % settings
clear data;

% average day per fish
data(:,1) = nanmean(compressibility(i_experiment_reps == er & i_group_tags == g,dn_hour == 1,1),2);
% average night per night
data(:,2) = nanmean(compressibility(i_experiment_reps == er & i_group_tags == g,dn_hour == 2,1),2);

shuffled = [nanmean(nanmean(compressibility(i_experiment_reps == er & ...
    i_group_tags == g,dn_hour == 1,2:end),2),3) ...
    nanmean(nanmean(compressibility(i_experiment_reps == er & ...
    i_group_tags == g,dn_hour == 2,2:end),2),3)];

data = [shuffled data];
data = data(:); 

% Grouping Variables
anova_group = [ones(size(data,1)/2,1) ; ones(size(data,1)/2,1)*2]; 
anova_time = [ones(size(data,1)/4,1) ; ones(size(data,1)/4,1)*2 ; ... 
    ones(size(data,1)/4,1) ; ones(size(data,1)/4,1)*2];
anova_development = ones(size(anova_group)); 
anova_experiment = repmat(i_experiment_tags(i_experiment_reps == er),4,1); 

% Comparison
[twa.p(:,1),~,twa.stats{er}] = anovan(data,...
    {anova_group,anova_time,anova_development,anova_experiment},...
    'display','off','model','full');


clear er anova_group anova_experiment anova_time t anova_development d data


%%  Deterministic Data 
scrap = [1 6 2 7 3 8 4 9 5 10];
scrap = repmat(scrap,1,step/length(scrap));
[~,~,compressibility_d] = compressSequenceNFast(scrap,sMax,nMax); 
compressibility_d = compressibility_d/step; 
compressibility_d = repmat(compressibility_d,length(i_group_tags(i_experiment_reps == er)),2); 

