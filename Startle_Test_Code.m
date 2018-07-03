%% Startle Test

%% Modules 

load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'i_experiment_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'time_window');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'lb'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'threads'); 
threads = threads(:,:,1); 

%% Modules Calc 

for f = 1:124
    set_token = i_experiment_tags(f); % used for each experiments sets settings
    counter = 1; 
    for bound = time_window{set_token}(1):time_window{set_token}(2)
        
        l_b = find(threads{f,2,1}(:,1) >= lb{set_token}(bound),1,'first');
        
        scrap(f,1:2,counter) = threads{f,1,1}(l_b-1:l_b,1); % modules 
        scrap(f,3:4,counter) = threads{f,2,1}(l_b-1:l_b,1) - lb{set_token}(bound); 
        counter = counter + 1; 
    end
end

%% Raw

load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'time_window');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'lb'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'i_experiment_tags');

load('D:\Behaviour\SleepWake\Re_Runs\Parameter_Extracted_data\Draft_1\170113_19_DATA.mat', 'delta_px_sq')
data{1} = delta_px_sq; 
load('D:\Behaviour\SleepWake\Re_Runs\Parameter_Extracted_data\Draft_1\170127_29_DATA.mat', 'delta_px_sq')
data{2} = delta_px_sq;
load('D:\Behaviour\SleepWake\Re_Runs\Parameter_Extracted_data\Draft_1\170210_09_DATA.mat', 'delta_px_sq')
data{3} = delta_px_sq;

for f = 1:3
    counter = 1;
    for bound = time_window{f}(1):time_window{f}(2)
        
        scrap{counter} = [scrap{counter} data{f}((lb{f}(bound) - ...
            (30*25) + 1):(lb{f}(bound) + (30*25)),:)];
        counter = counter + 1;
    end
end

plot([750 750],[0 70],'k','linewidth',3); 

