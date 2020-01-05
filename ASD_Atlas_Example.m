%% Atlas Example Figures 

num = 90; 
dot = 180; 

%% General Figure Formatting 
clf; hold on;
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); set(gca,'FontName','Calibri'); % Set Font
set(gca,'Xtick',[]); 
set(gca,'Ytick',[]); 
xlabel('Behaviour','Fontsize',32); 
ylabel('Brain','Fontsize',32); 
axis([1 num 1 num]); 

%% No Subtypes - V1 
scatter(1:num,repmat(num/2,num,1),...
    'k','filled','markerfacealpha',0.5)

%% No Subtypes - V2 
scatter(datasample(1:num,num,'Replace',false),...
    datasample(1:num,num,'Replace',false),...
    dot,'filled','k','markerfacealpha',0.5) 

%% All Subtypes (in 2 subtypes)
cmap = lbmap(2,'RedBlue');
scatter(normrnd(num/4,7,[num/2 1]),normrnd(num/4,7,[num/2 1]),dot,...
    'markerfacecolor',cmap(1,:), 'markeredgecolor', cmap(1,:),'markerfacealpha',0.5); 
scatter(normrnd(num/1.4,7,[num/2 1]),normrnd(num/1.4,7,[num/2 1]),dot,...
        'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5); 

%% Middle (2 subtypes + single points) 
clf; hold on; 
cmap = [lbmap(2,'RedBlue') ; [0 0 0]];
data(:,1) = [datasample(1:num,num/3,'Replace',false)' ; normrnd(num/4,3,[num/3 1]) ; ...    
normrnd(num/1.4,3,[num/3 1])];
data(:,2) = [datasample(1:num,num/3,'Replace',false)' ; normrnd(num/4,3,[num/3 1]) ; ...    
normrnd(num/1.4,3,[num/3 1])];

Mdl = fitgmdist(data,3); 
idx = cluster(Mdl,data); 

for i = 1:max(idx) 
    scatter(data(idx == i,1),data(idx == i,2),dot,...
        'markerfacecolor',cmap(i,:), 'markeredgecolor', cmap(i,:),'markerfacealpha',0.5); 
end 

%% Matrix Figure 
clf; hold on;
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); set(gca,'FontName','Calibri'); % Set Font
imagesc(data);
set(gca,'Xtick',[250 750]); 
set(gca,'TickLength',[0 0]);
set(gca,'XtickLabels',{'Behaviour','Brains'},'Fontsize',32); 
set(gca,'Ytick',[]); 
xlabel('Features','Fontsize',32); 
axis tight
cmap = flip(lbmap(9,'RedBlue'));
colormap(cmap); 

%% ASD 
data = [reshape((datasample(1:num,(num/3)*1000)),num/3,[]); ...
    [normrnd(num/4,10,[num/3 500]) [normrnd(num/2,10,[num/3 500])]] ...
    ; normrnd(num/1.4,10,[num/3 1000])];  
ylabel('Autism-associated Genes','Fontsize',32); 

%% Pharmacological 

data = rand(1000);

for i = 1:100:length(data)
    data(i:i+99,i:i+99) = data(i:i+99,i:i+99)*4; 
end 

ylabel('Pharmacological Compounds','Fontsize',32); 

%% Genetic 
sep = [sort(datasample(1:1000,10)) 1000]; 

data = rand(1000);
for i = 1:length(sep)-1
    data(sep(i):sep(i+1),sep(i):sep(i+1)) = data(sep(i):sep(i+1),sep(i):sep(i+1))+1; 
end 

ylabel('Genes','Fontsize',32); 
