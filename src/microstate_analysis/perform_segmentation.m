function [EEGp] = perform_segmentation(EEG, ALLEEG, options)
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
%                      Reject peaks X std over threshold
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
    EEG struct;
    ALLEEG struct = struct([]);
    options.remove_ocular_elecs logical = 1;
    options.select_Nmicro logical = 0;
    options.MinPeakDist = 10;
    options.Npeaks = 5000;
    options.GFPthresh = 1;
    options.Nmicrostates = 2:10;
    options.Nrepetitions = 50;
    options.algorithm = 'modkmeans';
    options.fitmeas = 'GEV'
    options.global_segmentation logical = 0;
end
if options.remove_ocular_elecs
    EEG.data = EEG.data(1:62,:);
    EEG.chanlocs = EEG.chanlocs(1:62);
end
% prepare the data for segmentation
if options.global_segmentation
    [~, ~, ~, EEGp] = pop_micro_selectdata(EEG,ALLEEG,'datatype', 'spontaneous', 'dataset_idx', [1:size(ALLEEG,2)],'MinPeakDist', options.MinPeakDist,'Npeaks', options.Npeaks,'GFPthresh', options.GFPthresh);
else
    EEGp = pop_micro_selectdata(EEG,[EEG],'datatype', 'spontaneous','MinPeakDist', options.MinPeakDist,'Npeaks', options.Npeaks,'GFPthresh', options.GFPthresh);
end
% perform the segmentation and create microstates prototypes
EEGp = pop_micro_segment(EEGp,'algorithm', options.algorithm,'Nmicrostates', options.Nmicrostates,'Nrepetitions', options.Nrepetitions, 'fitmeas', options.fitmeas);
% (selection of prototypes with eigenvectors in function 'segment' called by modkmeans)
% create a plot that helps you select the optimum number of microstates
if options.select_Nmicro
    EEGp = pop_micro_selectNmicro(EEGp, 'Measures', 'ALL');
end
end