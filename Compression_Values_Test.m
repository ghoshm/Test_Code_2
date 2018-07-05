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
    
    chunks{f,1}(end) = []; % remove the spare edge
    
    disp(num2str(f));
    
end
toc

clear f tc t

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

%%  Deterministic 
scrap = [1 6 2 7 3 8 4 9 5 10];
scrap = repmat(scrap,1,step/length(scrap));
[~,~,totSavings] = compressSequenceNFast(scrap,sMax,nMax); 

%% Random Sampling 

for f = 1:500
    a = datasample(1:2,1); 
    b = setdiff(1:2,a); 
    scrap = nan(1,step,'single');
    scrap(a:2:end) = datasample(1:5,step/2,'Replace',true);
    scrap(b:2:end) = datasample(6:10,step/2,'Replace',true);
    [~,~,totSavings(f)] = compressSequenceNFast(scrap,sMax,nMax);
end
