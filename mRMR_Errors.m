%% Collecting Data 

for t = 1:length(mRMR_ms) % for each comparison 
    scrap(t,1) = Mdl_loss{er,1}(t,mRMR_ms(er,t))*100; % collect percentage error
    scrap(t,2) = Mdl_loss{er,2}(t,mRMR_ms(er,t))*100; % collect percentage std 
end 

%% Null Classifiers
%http://blog.revolutionanalytics.com/2016/03/classification-models.html

n = size(mRMR_tw{er,1},1); % number of samples 
p = sum(mRMR_tw{er,1} == 1)/n; % prob group 

% Weighted Guess 
% null_w(1,1) = (1 - (p^2 + (1-p)^2)); % percentage error  
% null_w(1,2) = sqrt((null_w(1,1)*(1-null_w(1,1)))/n); % percentage std 
% null_w = null_w*100; 

% Majority Class 
null_m(1,1) = p*100;
null_m(1,2) = sqrt((p*(1-p))/n)*100; 

%% 
histcounts(mRMR_tw{er,1})

n = 28; 
p = 9/n; 

% Majority Class 
null_m(1,1) = p*100;
null_m(1,2) = sqrt((p*(1-p))/n)*100; 

%% Figure 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\mRMR_Comparisons\table.mat'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'cmap'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'cmap_2');
cmap{4} = [cmap{4}(2,:) ; cmap{4}(3,:) ; [1 0.5 0]] ;

% colormap 
n = 24; % number of colours 
CT = [linspace(cmap_2{1,1}(1,1),cmap_2{1,1}(2,1),n)'...
    linspace(cmap_2{1,1}(1,2),cmap_2{1,1}(2,2),n)'...
    linspace(cmap_2{1,1}(1,3),cmap_2{1,1}(2,3),n)']; 

cmap = [cmap{1} ; cmap{1} ; cmap{1} ; CT ; cmap{1} ; cmap{1};...  
    cmap{4}; cmap{4} ; cmap{4} ; cmap{6} ; cmap{7}]; 

% data 
clf; hold on; 
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
errorbar(data(:,3),data(:,4),'o','color',([1 1 1]*(1-(1/(9)^.5))),'linewidth',3); 
for i = 1:size(data,1)
    
    errorbar(i,data(i,1),data(i,2),'color',cmap(i,:),...
        'marker','o','linewidth',3);
end

% lines 
plot([1 29],[55 55],'color',cmap_2{1}(2,:),'linewidth',3); 
plot([30 38],[55 55],'color',[1 0.5 0],'linewidth',3); 
plot([39 45],[55 55],'color',cmap(39,:),'linewidth',3); 
plot([46 49],[55 55],'color',cmap(end,:),'linewidth',3); 

% labels 
text(12,57,'Wild Type','FontSize',32,'FontName','Calibri')
text(32,57,'HcrtR','FontSize',32,'FontName','Calibri')
text(39,57,'Melatonin','FontSize',32,'FontName','Calibri')
text(46.5,57,'PTZ','FontSize',32,'FontName','Calibri')

axis([0 size(data,1)+1 0 60]); 
set(gca,'XTick',[]);
xlabel('Classifier','Fontsize',32); 
ylabel('Classification Error (%)','Fontsize',32);