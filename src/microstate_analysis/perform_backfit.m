function [EEGp] = perform_backfit(EEGp, options)
%% Pipeline to perform microstate analysis
%
% **Usage:**
%   - perform_microstate_analysis(EEG)
%   - perform_microstate_analysis(..., Nmicrostates, [2:15])
%
% Input(s):
%    EEG = EEG struct
%
% Optional input parameter(s):
%    remove_ocular_elecs (default : 1) = 1 truncate electrode 63 and 64, 0 otherwise
%    select_Nmicro (default : 0) = 1 create the plot to select the number 
%                      of microstates,0 otherwise
%    MinPeakDist (default : 10) = input for pop_micro_selectdata - 
%                      Minimum Distance between GFP peaks in ms
%    Npeaks (default : 1000) = input for pop_micro_selectdata - 
%                      Number of GFP peaks per subject that enter the
%                      segmentation. Note that the maximum number of peaks
%                      is restricted to the minimum number of GFP peaks
%                      across subjects.
%   GFPthresh (default : 0) = input for pop_micro_selectdata - 
%                      Reject peaks over threshold
%   Nmicrostates (default : 2:10) = input for pop_micro_segment - 
%                      Number of microstates prototypes to consider
%   Nrepetitions (default : 50) = input for pop_micro_segment - 
%                      Number of times to repeat segmentation with random
%                      initialization
%   algorithm (default : modkmeans) = input for pop_micro_segment - see function
%   polarity (default : 0) = input for pop_micro_fit - 
%                      1 take polarity into account, 0 otherwise
%   perform_smoothing (default : 1) = perform smoothing where a minimum 
%                      duration is set for a microstate segment
%   minTime_smoothing (default : 30) input for pop_micro_smooth 
%                      only set if perform smoothing is 1 -
%                      Redristibute segments smaller than minTime (in ms) 
%                      to the next best fitting microstate. 
%   epoch (default : empty array) = input for pop_micro_stats -
%                      Timewindow of analysis (vector of timeframe indices). 
%                      If empty, all time samples will be used.
%
%
% Output(s):
%   EEGp = EEG struct with added microstate field containing segmentation
%         results and the stats
%
% Ruchella Kock, Leiden University, 22/08/2023 
%%
arguments
    EEGp struct;
    options.polarity logical = 0;
    options.perform_smoothing logical = 1;
    options.minTime_smoothing = 30;
    options.epoch = [];
    options.global_segmentation logical = 0;
end

% backfit the prototypes onto the data
EEGp = pop_micro_fit(EEGp, 'polarity', options.polarity);

% perform smoothing where a minimum duration is set for a microstate segment
if options.perform_smoothing
    % smooth the data otherwise the durations are too low
    EEGp = pop_micro_smooth(EEGp,'label_type', 'backfit','smooth_type', 'reject segments','minTime', options.minTime_smoothing,'polarity', options.polarity);
end

% calculate the stats over all_data
EEGp = pop_micro_stats(EEGp, 'label_type', 'backfit', 'epoch', options.epoch);
end