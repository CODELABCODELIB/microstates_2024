function [ALLEEG] = main_microstate_individual_level(ALLEEG,options)
%% Main function to run the microstate analysis at the individual level
%  this function runs it for movie data, resting state data or sequences
%  where there are tap data
%
% **Usage:**
%   - [ALLEEG] = main_microstate_individual_level(ALLEEG)
%               - main_microstate_individual_level(...,'start_run', 20)
%               - main_microstate_individual_level(...,'end_run', 25)
%               - main_microstate_individual_level(...,'plot_data', 0)
%
% Input(s):
%    ALLEEG = all EEG structs for all participants without microstate field
%
% Optional input parameter(s):
%   tap_only (logical, default : 1) = 1 microstate analysis on tap data
%   movie_data (logical, default : 0) = 1 microstate analysis on movie data
%   resting_state (logical, default : 0) = 1 microstate analysis on resting
%           state data at the individual level
%   start_run (default : 1) = index to start with in ALLEEG
%   end_run (default : size(ALLEEG,1)) = index to end with in ALLEEG
%   plot_data (logical, default : 0) = create plot of the microstates for
%         each participant
%
% Output(s):
%    ALLEEG = all EEG structs for all participants with microstate field
%       if 'tap_only' is selected then ALLEEG(:,3) includes the indexes
%       that were used where tap data is present
%
% Requires:
%   prepare_EEG_w_taps_only.m
%   perform_segmentation.m
%   perform_backfit.m
%   find_movie_passive_event.m
%   perform_microstate_analysis.m 
%   prepare_EEG_w_resting_state.m
%   plot_microstates_topoplots.m
%
% Ruchella Kock, Leiden Unviersity, 12/09/2023
arguments
    ALLEEG
    options.plot_data = 0;
    options.tap_only = 1;
    options.movie_data = 0;
    options.resting_state = 0;
    options.start_run = 1;
    options.end_run = size(ALLEEG,1)
    options.file = '';
    options.save_results logical = 1;
    options.save_path = '.';
    options.parameters = [];
    options.threshold = 0; 
end
% run for every participant
for pp=options.start_run:options.end_run
    EEG = ALLEEG{pp,2};
    if options.tap_only
        [EEGtap, indexes] = prepare_EEG_w_taps_only(EEG);
        [EEGtap] = perform_segmentation(EEGtap);
        [EEGtap] = perform_backfit(EEGtap);
        %         [EEGtap] = perform_backfit(EEGtap,'epoch', [indexes{:}]);
        ALLEEG{pp,2} = EEGtap;
        ALLEEG{pp,3} = indexes;
    end
    if options.movie_data
        [~,~,movie_present,latencies] = find_movie_passive_event(EEG);
        [EEGmoviestats] = perform_microstate_analysis(EEG,'epoch', latencies(1):latencies(end));
        ALLEEG{pp,4} = EEGmoviestats;
    end
    if options.resting_state
        [EEG_EO,EEG_EC] = prepare_EEG_w_resting_state(EEG);
        % eyes open
        [EEG_EO] = perform_segmentation(EEG_EO);
%         [EEG_EO] = perform_backfit(EEG_EO);
        % eyes close
        [EEG_EC] = perform_segmentation(EEG_EC);
%         [EEG_EC] = perform_backfit(EEG_EC);
        ALLEEG{pp,5} = EEG_EO;
        ALLEEG{pp,6} = EEG_EC;
    end
    % plot the maps
    if options.plot_data
        plot_microstates_topoplots(EEG);
    end
end
%% save cdresults
if options.save_results
    save(sprintf('%s/%s',options.save_path,options.file), 'ALLEEG', '-v7.3')
end