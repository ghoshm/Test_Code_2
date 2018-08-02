% Test 

%% mRMR Ideas 

% Clear data before processing 
clear comps* er set_token mRMR* Mdl*
comps = 1000; 
er = 1; 
set_token = find(experiment_reps == er,1,'first'); % settings

% Data to Save 
clearvars -except comps comps_v er set_token Mdl_loss mRMR_data mRMR_tw mRMR_ms mRMR_tsne 

%% Specify which data you want 
% Development 

% WT Keep Days 
mRMR_data{er,1}(249:end,:) = []; 
mRMR_tw{er,1}(249:end,:) = []; 
mRMR_tw{er,1}(125:end) = 2; 

% WT Keep Nights 
mRMR_data{er,1}(1:248,:) = []; 
mRMR_tw{er,1}(249:end,:) = []; 
mRMR_tw{er,1}(125:end) = 2; 

% REMOVE NIGHTS FROM OTHERS 
mRMR_data{er,1}(2:2:end,:) = []; 
mRMR_tw{er,1}(2:2:end,:) = []; 

% REMOVE DAYS FROM OTHERS 
mRMR_data{er,1}(1:2:end,:) = []; 
mRMR_tw{er,1}(1:2:end,:) = []; 

% Evening/Morning 
mRMR_tw{1,1}(mRMR_tw{1,1} <= 7,:) = 1;
mRMR_tw{1,1}(mRMR_tw{1,1} > 7,:) = 2; 

% Early vs Late Night 
mRMR_data{er,1}(mRMR_tw{er,1} <= 14,:) = []; % keep only some data
mRMR_tw{er,1}(mRMR_tw{er,1} <= 14,:) = []; % keep only some data
mRMR_tw{1,1}(mRMR_tw{1,1} <= 19,:) = 1;
mRMR_tw{1,1}(mRMR_tw{1,1} > 19,:) = 2; 

% Hours Averaging  

% Startles 
mRMR_tw{1,1}(mRMR_tw{1,1} == 15,:) = 1; 
mRMR_tw{1,1}(mRMR_tw{1,1} ~= 1,:) = 2;

%% Corrupt WT Labels
    % To test Hcrt Data  

tic 
counter = 1; % start a counter (counts comparisons)
for r = 1:10 % for each repeat 

    % Corrupt the labels
    mRMR_tw{er,1} = repmat(randsrc(length(mRMR_tw{er,1})/2,1,[1 2 3 ; ...
        histcounts(i_group_tags(i_experiment_reps == 2))/sum(i_experiment_reps == 2)]),2,1);

    %  % Pairwise Comparisons
    for g_one = min(mRMR_tw{er,1}):max(mRMR_tw{er,1}) % for each group
        for g_two = (g_one + 1):max(mRMR_tw{er,1}) % for each comparison

            % mRMR
            [comps_v{er,1}(counter,:)] = mrmr_miq_d(...
                zscore(mRMR_data{er,1}(mRMR_tw{er,1} == g_one | mRMR_tw{er,1} == g_two,:)),...
                mRMR_tw{er,1}(mRMR_tw{er,1} == g_one | mRMR_tw{er,1} == g_two,:),comps);

            % Classifiers
            for s = 1:comps % for each comp sequence
                % Fit a linear classifier as you add features
                % Using 10 fold cross validation
                % Hold 10% of mRMR_data back by default
                Mdl = fitcdiscr(...
                    zscore(mRMR_data{er,1}(mRMR_tw{er,1} == g_one | mRMR_tw{er,1} == g_two,...
                    comps_v{er,1}(counter,1:s))),...
                    mRMR_tw{er,1}(mRMR_tw{er,1} == g_one | mRMR_tw{er,1} == g_two,:),...
                    'DiscrimType','linear','CrossVal','on');
                Mdl_loss{er,1}(counter,s) = kfoldLoss(Mdl);
                Mdl_loss{er,2}(counter,s) = nanstd(kfoldLoss(Mdl,'Mode','individual'));
                %disp(num2str(s));
            end

            % Minimal Feature Space
            % Note 180509: could also use islocalmin
            if er == 1 % for the WT data
                mRMR_ms(er,counter) = find(smooth(Mdl_loss{er,1}(counter,:),3) == ...
                    min(smooth(Mdl_loss{er,1}(counter,:),3)),1,'first');
            else
                %try
                %mRMR_ms(er,counter) = find(Mdl_loss{er,1}(counter,:) < 0.05,1,'first');
                %catch
                mRMR_ms(er,counter) = find(smooth(Mdl_loss{er,1}(counter,:),3) == ...
                    min(smooth(Mdl_loss{er,1}(counter,:),3)),1,'first');
                %end
            end

            counter = counter + 1; % add to counter (counts comparisons)
        end
    end
    disp(horzcat('Finished mRMR Comparisons ',num2str(er),' of ',...
        num2str(max(experiment_reps)))); % report progress

end
toc 

%% NEWER MESSING AROUND 
data = gCount_norm{1,1}(comps_v{er,1}(2,1:mRMR_ms(er,2)),[4 6],i_experiment_reps == er);

data = reshape(data,mRMR_ms(er,2),[]); 

% Remove Hets 
data(:,mRMR_tw{er,1} == 2) = []; 
mRMR_tw{er,1}(mRMR_tw{er,1} == 2) = []; 

% Separting day/night 
mRMR_tw{er,1}(2:2:end,:) = mRMR_tw{er,1}(2:2:end,:)*2; 
mRMR_tw{er,1}(mRMR_tw{er,1} == 6) = 4; 

% scatter 
gscatter(data(1,:),data(2,:),mRMR_tw{er,1});

% line plot 
temp = []; 
for s = 1:size(data,1) 
    temp = [temp grpstats(data(s,:),mRMR_tw{er,1},'mean')]; 
end 
[~,O] = sort(temp(1,:) - temp(2,:));

clf; hold on; 
plot(temp(1,O)); 
plot(temp(2,O)); 

% Motifs Involved 
idx = temp(1,:) > temp(2,:); 
temp = grammar_mat{1,1}(comps_v{er,1}(2,idx==0),:); 
temp = temp(:); 
temp(isnan(temp)) = [];
clf; hold on; 
plot(histcounts(temp,'normalization','pdf'));

temp = grammar_mat{1,1}(comps_v{er,1}(2,idx==1),:); 
temp = temp(:); 
temp(isnan(temp)) = [];
plot(histcounts(temp,'normalization','pdf'));

% Motif Lengths 
temp = grammar_mat{1,1}(comps_v{er,1}(2,idx==0),:); 
clf; hold on; 
histogram(sum(isnan(temp)==0,2),'Normalization','probability'); 

temp = grammar_mat{1,1}(comps_v{er,1}(2,idx==1),:); 
histogram(sum(isnan(temp)==0,2),'Normalization','probability'); 

% pca
[coeff,score,~,~,explained,~] = pca(data'); 
clf; 
gscatter(score(:,1),score(:,2),mRMR_tw{er,1}); 

% tsne 
mRMR_tsne{er,1} = tsne(data',...
    'Algorithm','exact','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',5,...
    'Perplexity',30,'Standardize',1,'Verbose',1);

clf; hold on;
for g = [1 3]
    x = double(mRMR_tsne{er,1}(mRMR_tw{er,1} == g,1)); 
    y = double(mRMR_tsne{er,1}(mRMR_tw{er,1} == g,2)); 
    k = boundary(x,y);  
    patch(x(k),y(k),cmap{set_token}(g,:),'EdgeColor',cmap{set_token}(g,:),'FaceAlpha',.5);
    scatter(mRMR_tsne{er,1}(mRMR_tw{er,1} == g,1),mRMR_tsne{er,1}(mRMR_tw{er,1} == g,2),...
        'markerfacecolor',cmap{set_token}(g,:),...
        'markeredgecolor',cmap{set_token}(g,:));
end 

gscatter(mRMR_tsne{er,1}(:,1),mRMR_tsne{er,1}(:,2),mRMR_tw{er,1}); 

%% "Leave one out" Classifiers 

comps = 250; % number of motifs to find 

% Pairwise Comparisons
tic
counter = 1;
for g_one = min(mRMR_tw{er,1}):max(mRMR_tw{er,1}) % for each group (hour)
    
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
disp(horzcat('Finished mRMR Comparisons ',num2str(er),' of ',...
    num2str(max(experiment_reps)))); % report progress
toc
    
%% mean inactive module length (frames)
ibl = grpstats(sleep_cells(:,3),idx_numComp_sorted{2,1},'mean');
ibl(1) = []; % remove NaN's
ibl(end) = 250; % crop longest module to 10s (hard coded)

%% tSNE Calculation 
er = 3; 
set_token = find(experiment_reps == er,1,'first'); % settings
    
motifs = [];
for c = 1:size(comps_v{er,1},1) % for each comparison
    motifs = [motifs comps_v{er,1}(c,1:mRMR_ms(er,c))];
end 
motifs = unique(motifs); % keep only unique motifs    

data = gCount_norm{1,1}(motifs,...
    1,...
i_experiment_reps == er); % take DAY data 

data = squeeze(data); % reshape to motifs x fish 

[coeff,score,~,~,explained,~] = pca(data'); % pca 

mRMR_tsne{er,1} = tsne(data',...
    'Algorithm','exact','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
    'Perplexity',30,'Standardize',1,'Verbose',1); % tSNE

%% Leave One Out Classifiers Figure 

figure;
clear legend_cols legend_cell

% Motifs
subplot(1,3,1); % subplot
hold on; set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
b = 1; % baseline counter (plots from bottom to top)
for s = flip(comps_v{er,1}(:,1))' % for each motif
    seq = grammar_mat{1,1}(s,:); % find this motifs module sequence
    seq(isnan(seq)) = []; % remove nan values
    a = 1; % start a counter (frames)

    for t = 1:length(seq) % for each module in the sequence
        if seq(t) <= numComp(1) % for the inactive modules
            plot([a (a+ibl(seq(t)))],[b b],...
                'color',cmap_cluster_merge(seq(t),:),'linewidth',5); % plot
            a = a + ibl(seq(t)); % add to time
        else % for the active modules
            plot(a:(a+length(nanmean(bouts{1,seq(t)-numComp(1)}))+1),...
                [b ((nanmean(bouts{1,seq(t)-numComp(1)})/28)+b) b],...
                'color',cmap_cluster_merge(seq(t),:),'linewidth',5); % plot
            a = a + length(nanmean(bouts{1,seq(t)-numComp(1)})) + 1; % add to time
        end
    end

    b = b + 1; % add to baseline counter
end

x_lims = xlim; 
x_lims = ceil(x_lims(2)/25)*25; 
axis([1 x_lims .5 b]); % hard coded axis
set(gca,'XTick',[1 x_lims/2 x_lims]); % hc x axis ticks
set(gca,'XTickLabels',[1/25 (x_lims/25)/2 x_lims/25]); % hc x axis labels
xlabel('Time (Seconds)','Fontsize',32);
ylabel('Comparison','Fontsize',32);
set(gca,'YTick',1:max(i_group_tags(i_experiment_reps == er)));
set(gca,'YTickLabels',flip(geno_list{set_token}.colheaders),'Fontsize',32); 

% Legend
% for s = 1:length(cmap_cluster_merge) % for each module
%     if s <= numComp(1) % for the inactive modules
%         scatter(max(xlim)-40,(s/5)+(size(comps_v{er,1},1)-.5),90,'markerfacecolor',cmap_cluster_merge(s,:),...
%             'markeredgecolor','k');
%     else % for the active modules
%         scatter(max(xlim)-10,((s-numComp(1))/5)+(size(comps_v{er,1},1)-.5),90,'markerfacecolor',cmap_cluster_merge(s,:),...
%             'markeredgecolor','k');
%     end
% end
% text(max(xlim)-95,(size(comps_v{er,1},1)+.75),'Module','FontName','Calibri','Fontsize',24);

% Scores 
subplot(1,3,2); 
hold on; set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);

data = gCount_norm{1,1}(comps_v{er,1}(:,1),1,i_experiment_reps == er); % DAY 1 
data = squeeze(data); % motifs x fish
data = flip(data); % flipped  
shift = abs(min(data(:))) + 1; % shift so that min(data(:)) == 1 (for Log axis) 
plot([shift shift],[.5 size(comps_v{er,1},1)+1],'-k','linewidth',1.5); 

data = data + shift; % shift data so that min(data(:)) == 1 

for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group 
    spread_cols = plotSpread(data(:,i_group_tags(i_experiment_reps == er) == g)',...
        'xyOri','flipped','spreadWidth',1/max(i_group_tags(i_experiment_reps == er)),...
        'distributionColors',...
        cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'XValues',...
        (1:size(data,1)) + (max(i_group_tags(i_experiment_reps == er))-1)/10 - (g-1)/10);
    set(findall(gca,'type','line'),'markersize',9); % change marker size 
    
    legend_cols(g,:) = errorbar(nanmean(data(:,i_group_tags(i_experiment_reps == er) == g),2),...
        ((1:size(data,1)) + (max(i_group_tags(i_experiment_reps == er))-1)/10 - (g-1)/10),...
        nanstd(data(:,i_group_tags(i_experiment_reps == er) == g)'),...
        'horizontal','o','linewidth',3,'color',cmap{set_token}(g,:),'capsize',9);
    legend_cell{g} = horzcat(' ',geno_list{set_token}.colheaders{g},', n = ',...
        num2str(sum(i_group_tags(i_experiment_reps == er) == g)));
    
end

% Legend
[~,icons,plots,~] = legend(legend_cols,legend_cell,'Location','best','FontSize',24);
legend('boxoff'); 
set(plots,'LineWidth',3); 

axis([min(data(:)) max(data(:)) .5 b]); % hard coded axis
set(gca,'XTick',[shift round(max(data(:)))]);
set(gca,'XTickLabels',string(num2cell(round([shift max(data(:))] - shift))));
%set(gca,'XScale','log'); 
xlabel('Z-Score','Fontsize',32);
set(gca,'YTick',[]);

% tSNE  
subplot(1,3,3); 
hold on; set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);

for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group 
    x = double(mRMR_tsne{er,1}(i_group_tags(i_experiment_reps == er) == g,1)); 
    y = double(mRMR_tsne{er,1}(i_group_tags(i_experiment_reps == er) == g,2)); 
    k = boundary(x,y);  
    patch(x(k),y(k),cmap{set_token}(g,:),'EdgeColor',cmap{set_token}(g,:),'FaceAlpha',.5);
end 

for g = 1:max(i_group_tags(i_experiment_reps == er))
    scatter(mRMR_tsne{er,1}(i_group_tags(i_experiment_reps == er) == g,1),...
        mRMR_tsne{er,1}(i_group_tags(i_experiment_reps == er) == g,2),60,...
        'markerfacecolor',cmap{set_token}(g,:),...
        'markeredgecolor',cmap{set_token}(g,:)); 
end 
set(gca,'XTick',[]); set(gca,'YTick',[]); 
xlabel('tSNE 1','Fontsize',32); 
ylabel('tSNE 2','Fontsize',32); 
axis tight; 

%% Hours Data 

data = mRMR_data{1,1}(:,comps_v{er,1}(:,1:mRMR_ms)); 

temp = [];
for s = 1:size(data,2) 
    temp = [temp grpstats(data(:,s),mRMR_tw{1,1},'mean')];
end 

mRMR_tsne{er,1} = tsne(data,...
    'Algorithm','barneshut','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
    'Perplexity',30,'Standardize',1,'Verbose',1); 
gscatter(mRMR_tsne{er,1}(:,1),mRMR_tsne{er,1}(:,2),mRMR_tw{er,1})

[coeff,score,~,~,explained,~] = pca(data); % pca
gscatter(score(:,1),score(:,2),mRMR_tw{er,1})

%% No mRMR? 

data = gCount_norm{1,1}(:,:,i_experiment_reps == 1);
data = nanmean(data,3);
[coeff,score,~,~,explained,~] = pca(data); % pca

mRMR_tsne{er,1} = tsne(data,...
    'Algorithm','barneshut','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
    'Perplexity',30,'Standardize',1,'Verbose',1); 

gscatter(mRMR_tsne{er,1}(:,1),mRMR_tsne{er,1}(:,2),mRMR_tw{er,1})
%% Learning about the "Languauge" 

% Modules Used
figure; hold on;
for c = 1:size(comps_v{er,1},1)
    data = grammar_mat{1,1}(comps_v{er,1}(c,1:mRMR_ms(er,c)),:);
    data = data(:);
    data(isnan(data)) = [];
    plot(histcounts(data,'BinMethod','integers','normalization','probability'),...
        'color',cmap{6}(c,:),'linewidth',3);
end

% Motif Lengths 
figure; hold on;
for c = 1:size(comps_v{er,1},1)
    data = grammar_mat{1,1}(comps_v{er,1}(c,1:mRMR_ms(er,c)),:); 
    data = sum(isnan(data) == 0,2); 
    plot(histcounts(data,'BinMethod','integers','normalization','probability'),...
        'color',cmap{6}(c,:),'linewidth',3);
end

%% Data to scatter & Motif order 
data = gCount_norm{1,1}(comps_v{er,1}(1:mRMR_ms(er,1)),...
    [1 2],...
i_experiment_reps == er); 

temp = squeeze(data(:,1,:))' - squeeze(data(:,2,:))'; 
[~,O] = sort(nanmean(temp)); % mean difference (ascending)   
temp = temp(:,O); 
data = data(O,:,:); % sort by mean difference 

clf; hold on; 
plotSpread(temp); % difference 

clf; hold on; 
errorbar(nanmean(data(:,1,:),3),nanstd(squeeze(data(:,1,:))')/sqrt(124));
errorbar(nanmean(data(:,2,:),3),nanstd(squeeze(data(:,2,:))')/sqrt(124),'r')

data = reshape(data,size(data,1),[])'; 
[coeff,score,~,~,explained,~] = pca(data); % pca

clf; hold on; 
scatter(score(1:124,1),score(1:124,2),'filled');
scatter(score(125:end,1),score(125:end,2),'filled');

mRMR_tsne{er,1} = tsne(data,...
    'Algorithm','exact','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
    'Perplexity',30,'Standardize',1,'Verbose',1);

clf; hold on; 
scatter(mRMR_tsne{er,1}(1:124,1),mRMR_tsne{er,1}(1:124,2),'filled');
scatter(mRMR_tsne{er,1}(125:end,1),mRMR_tsne{er,1}(125:end,2),'filled');

%% Embedding Melatonin in the WT Space 
er = 3; 
set_token = find(experiment_reps == er,1,'first'); % settings
    
data = mRMR_data{er,1}(:,comps_v{1,1}(1,1:mRMR_ms(1)));
data = data(:,O); 

clf; hold on;
col = 1;
for g = 1:max(mRMR_tw{er,1})
    temp = data(mRMR_tw{er,1} == g,:);
    
    scatter(temp(1:2:end,1),temp(1:2:end,2),'filled',...
        'markerfacecolor',cmap_2{set_token}(col,:));
    scatter(temp(2:2:end,1),temp(2:2:end,2),'filled',...
        'markerfacecolor',cmap_2{set_token}(col+1,:));
    
    col = col + 2;
    
end

[coeff,score,~,~,explained,~] = pca(data); % pca

clf; hold on;
col = 1;
for g = 1:max(mRMR_tw{er,1})
    temp = score(mRMR_tw{er,1} == g,:);
    
    scatter(temp(1:2:end,1),temp(1:2:end,2),'filled',...
        'markerfacecolor',cmap_2{set_token}(col,:));
     scatter(temp(2:2:end,1),temp(2:2:end,2),'filled',...
        'markerfacecolor',cmap_2{set_token}(col+1,:));
    
    col = col + 2;
    
end

scrap = tsne(data,...
    'Algorithm','exact','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
    'Perplexity',30,'Standardize',1,'Verbose',1);

clf; hold on;
col = 1;
for g = 1:max(mRMR_tw{er,1})
    temp = scrap(mRMR_tw{er,1} == g,:);
    
    scatter(temp(1:2:end,1),temp(1:2:end,2),'filled',...
        'markerfacecolor',cmap_2{set_token}(col,:),...
        'markerfacealpha',0.5);
    plot(nanmean(temp(1:2:end,1)),nanmean(temp(1:2:end,2)),...
        'x','color',cmap_2{set_token}(col,:),'linewidth',3,'MarkerSize',15)
    
    scatter(temp(2:2:end,1),temp(2:2:end,2),'filled',...
        'markerfacecolor',cmap_2{set_token}(col+1,:),...
        'markerfacealpha',0.5);
    plot(nanmean(temp(2:2:end,1)),nanmean(temp(2:2:end,2)),...
        'x','color',cmap_2{set_token}(col+1,:),'linewidth',3,'MarkerSize',15)
    col = col + 2;
end

%% Plotting the Melatonin Space 

er = 3; 
set_token = find(experiment_reps == er,1,'first'); % settings
    
motifs = [];
for c = 1:size(comps_v{er,1},1) % for each comparison
    motifs = [motifs comps_v{er,1}(c,1:mRMR_ms(er,c))];
end 
motifs = unique(motifs);   

data = gCount_norm{1,1}(motifs,...
    1,...
i_experiment_reps == er); % take day data 

data = squeeze(data); % reshape to motifs x fish 

for m = 1:size(data,1) % for each motif 
    motif_sorting(m,:) = grpstats(data(m,:),...
        i_group_tags(i_experiment_reps == er),'mean');
end 
[~,idx] = max(motif_sorting'); 
[~,O] = sort(idx);
motif_sorting = motif_sorting(O,:);
idx = idx(O); 
data = data(O,:); 

for g = 1:max(i_group_tags(i_experiment_reps == er)) 
    [~,O] = sort(motif_sorting(idx == g,g),'descend'); 
    temp = data(idx == g,:); 
    temp = temp(O,:); 
    data(idx == g,:) = temp; 
end

clf; hold on;
for c = 1:max(idx)
    
    for g = 1:max(i_group_tags(i_experiment_reps == er))
        %
        %             plotSpread(data(idx == c,i_group_tags(i_experiment_reps == er) == g)',...
        %                 'distributionColor',cmap{set_token}(g,:),'xyOri','flipped')
        %         %     set(findall(gca,'type','line'),'markersize',9); % change
        
        plot(nanmean(data(idx ==c,i_group_tags(i_experiment_reps == er) == g),2),find(idx == c),...
            'color',cmap{set_token}(g,:),'linewidth',1.5);
        %             errorbar(nanmean(data(idx == c,i_group_tags(i_experiment_reps == er) == g),2),...
        %                 find(idx == c),nanstd(data(idx == c,i_group_tags(i_experiment_reps == er) == g)')/...
        %                 sum(i_group_tags(i_experiment_reps == er) == g),'horizontal','.','color',...
        %                 cmap{set_token}(g,:),'linewidth',3);
    end
    
    plot(nanmean(data(idx ==c,i_group_tags(i_experiment_reps == er) == c),2),find(idx == c),...
        'color',cmap{set_token}(c,:),'linewidth',5);
    
end

mRMR_tsne{er,1} = tsne(data',...
    'Algorithm','exact','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
    'Perplexity',30,'Standardize',1,'Verbose',1);

clf; hold on; 
for g = 1:max(i_group_tags(i_experiment_reps == er))
    scatter(mRMR_tsne{er,1}(i_group_tags(i_experiment_reps == er) == g,1),...
        mRMR_tsne{er,1}(i_group_tags(i_experiment_reps == er) == g,2),60,...
        'markerfacecolor',cmap{set_token}(g,:),...
        'markeredgecolor',cmap{set_token}(g,:)); 
end 

clf; hold on; 
for g = 1:max(i_group_tags(i_experiment_reps == er))
    x = double(mRMR_tsne{er,1}(i_group_tags(i_experiment_reps == er) == g,1)); 
    y = double(mRMR_tsne{er,1}(i_group_tags(i_experiment_reps == er) == g,2)); 
    k = boundary(x,y);  
    patch(x(k),y(k),cmap{set_token}(g,:),'EdgeColor',cmap{set_token}(g,:),'FaceAlpha',.5);
end 

[coeff,score,~,~,explained,~] = pca(data'); % pca
clf; hold on; 
for g = 1:max(i_group_tags(i_experiment_reps == er))
    scatter(score(i_group_tags(i_experiment_reps == er) == g,1),...
       score(i_group_tags(i_experiment_reps == er) == g,2),90,...
        'markerfacecolor',cmap{set_token}(g,:),...
        'markeredgecolor',cmap{set_token}(g,:)); 
end 


%% PTZ Data 
er = 4; 
set_token = find(experiment_reps == er,1,'first'); % settings

data = gCount_norm{1,1}(comps_v{er,1}(1:mRMR_ms(er,1)),...
    [1 2],...
i_experiment_reps == er); 

clf; hold on;
for g = 1:max(i_group_tags(i_experiment_reps == er))
    
    errorbar(nanmean(data(:,1,i_group_tags(i_experiment_reps == er) == g),3),...
        nanstd(squeeze(data(:,1,i_group_tags(i_experiment_reps == er) == g))'),'color',...
        cmap{set_token}(g,:));
    
end
