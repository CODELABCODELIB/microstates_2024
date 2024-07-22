function [h_values,p_values,test_values,stats_all] = main_statistical_analysis_prob(ALLEEG,options)
%% Main function to run the stat testing of the jid-states for all participants
%
% **Usage:**
%   - [h_values,p_values,test_values,stats_all] = main_statistical_analysis(ALLEEG)
%                - main_statistical_analysis(...,'start_run', 20)
%                - main_statistical_analysis(...,'end_run', 25)
%                - main_statistical_analysis(...,'save_results', 0)
%                - main_statistical_analysis(...,'test', 'ttest')
%
% Input(s):
%    ALLEEG struct = all EEG structs for all participants with microstate field
%           (shape : number of participants x 2) where the first column
%           contains the participant number (unused in function) and second
%           the EEG struct for each participant WITH microstates field
%
% Optional input parameter(s):
%   save_results logical (default : 1) = 1 save each jid-state and boot
%           jid-states, save the h,p and stats, 0 otherwise
%   save_path (default : current folder) = path to save results
%   verbose logical (default : 1) = 1 show which state and participant
%           is being processed, 0 otherwise
%   start_run (default : 1) = index to start with in ALLEEG
%   end_run (default : size(ALLEEG,1)) = index to end with in ALLEEG
%   test (string, default : 'signrank') = chosen test options include
%         'ttest' or 'signrank'
%
% Output(s):
%   h_values = result of the hypothesis test
%   p_values = p-value of the test
%   test_values = t-value or sign value
%   stats_all (cell array) =  Test statistics
%   (shape of each output: 2500 x max number of states (11) x number of participants)
%
% Requires:
%   get_microstate_labels.m
%   create_jid_state.m
%   create_jid_state_boots.m
%   get_test_statistic.m
%   create_boot_test_val.m
%   MCC_jids.m
%
% Ruchella Kock, Leiden Unviersity, 12/09/2023
arguments
    ALLEEG;
    options.test = 'prctile';
    options.save_results logical = 1;
    options.save_path = '.';
    options.file = '';
    options.verbose logical = 1;
    options.start_run = 1;
    options.end_run = size(ALLEEG,1);
    options.parameter = 1;
    options.pool_method = 'mean';
    options.threshold = 0; 
end
%% initialize arrays
h_values = zeros(2500,11,size(ALLEEG,1));
p_values = zeros(2500,11,size(ALLEEG,1));
test_values = zeros(2500,11,size(ALLEEG,1));
stats_all = cell(2500,11,size(ALLEEG,1));
%%
tic
% repeat the analysis per participant
for pp=options.start_run:options.end_run
    EEG = ALLEEG{pp,2};
    [taps] = find_taps(EEG, ALLEEG{pp,3}); % tap indexes
    [dt_dt,~] = calculate_ITI_K_ITI_K1(taps, 'shuffle', 0); %
    microstates = EEG.microstate.fit.labels; % all labels (size equal to EEG.data)
    n_states = max(microstates);
    
    [jid_state_original,probabilities] = original_state_prob(taps,microstates,dt_dt,n_states, 'duration', options.parameter, 'pool_method', options.pool_method);
    [jid_state_boot,probabilities_boot] = boot_state_prob(taps,microstates,dt_dt,n_states, 'duration', options.parameter, 'pool_method', options.pool_method);
    
    for state=1:n_states
        if options.verbose
            fprintf('Calculating state %d - Sub %d \n',state, pp)
        end
        parfor bin=1:2500
            [h_values(bin,state,pp),p_values(bin,state,pp),stats] = get_test_statistic(jid_state_boot{state}(bin,:),jid_state_original{state}(bin),'test', options.test);
            if strcmp(options.test, 'ttest')
                test_values(bin,state,pp) = stats.tstat;
            elseif strcmp(options.test, 'signrank')
                test_values(bin,state,pp) = stats.sign;
            else
                test_values(bin,state,pp) = stats;
            end
            stats_all{bin,state,pp} = stats;
        end
%         if sum(test_values(:,state,pp))
%             [boot_data] = create_boot_test_val(test_values(:,state,pp));
%             [cluster_states{state}] = MCC_jids(test_values(:,state,pp),boot_data);
%         end
    end
    if options.save_results
        save(sprintf('%s/%s_jid_state_boots_%d.mat',options.save_path, options.file, pp), 'jid_state_boot')
        save(sprintf('%s/%s_jid_state_%d.mat',options.save_path,  options.file,pp), 'jid_state_original')
%         save(sprintf('%s/%s_clustered_res_%d.mat',options.save_path, options.file,pp), 'cluster_states')
        save(sprintf('%s/%s_probabilities_%d.mat',options.save_path, options.file,pp), 'probabilities')
        save(sprintf('%s/%s_probabilities_boot_%d.mat',options.save_path,options.file, pp), 'probabilities_boot')
        
    end
end
if options.save_results
    save(sprintf('%s/h_values_%s_%d_to_%d.mat',options.save_path, options.file,options.start_run,options.end_run), 'h_values')
    save(sprintf('%s/p_values_%s_%d_to_%d.mat',options.save_path,options.file,options.start_run,options.end_run), 'p_values')
    save(sprintf('%s/test_values_%s_%d_to_%d.mat',options.save_path,options.file,options.start_run,options.end_run), 'test_values')
end
toc
end