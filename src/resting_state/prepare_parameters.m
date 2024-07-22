function [params_per_pps] = prepare_parameters(data_path,type,search_condition)
%% Perform regression between JID 1 month before measurement and microstate parameter
%
% **Usage:** [params_per_pps] =  prepare_parameters(data_path,type,search_condition)
%        - rest_correlation(..., search_condition,'EEG2')
%
%  Input(s):
%   - data_path = data_path string path to the backfit files (day 1 and day 2)
%   - type = string 'EO' Eyes open or 'EC' Eyes closed
%
%  Optional Input(s):
%   - search_condition (Default:'EEG') = String or char array file name
%       common accross all subjects containing the microstate stats (e.g. for
%       file EEG_EC_1.mat then 'EEG' is the search condition)
%
%  Output(s):
%   - params_per_pps = cell with microstates parameters for each subject 
%       params_per_pps{pp,1} = Duration;
%       params_per_pps{pp,2} = Occurence;
%       params_per_pps{pp,3} = Coverage;
%       params_per_pps{pp,4} = prototypes;
%       params_per_pps{pp,5} = subject folder name;
%       params_per_pps{pp,6} = GEV;
%
% Author: R.M.D. Kock
arguments
    data_path;
    type;
    search_condition = 'EEG2';
end
folder_contents = dir(data_path);
files = {folder_contents.name};
all_files = files(cellfun(@(x) contains(x, sprintf('%s_%s',search_condition,type)),files));
params_per_pps = cell(length(all_files),3);
for pp=1:length(all_files)
    % file selection
    load(all_files{pp})
    %     params_per_pps(pp,1) = {EEG_EO.microstate.stats.Duration(1:end-1)};
    params_per_pps(pp,1) = {EEG_EO.microstate.stats.Duration(1:end)};
    params_per_pps(pp,2) = {EEG_EO.microstate.stats.Occurence(1:end)};
    params_per_pps(pp,3) = {EEG_EO.microstate.stats.Coverage(1:end)};
    params_per_pps(pp,4) = {EEG_EO.microstate.prototypes(:,1:end)};
    params_per_pps(pp,5) = {EEG_EO.filepath(end-3:end)};
    params_per_pps(pp,6) = {EEG_EO.microstate.stats.GEV(1:end)};
    clearvars EEG_EO;
end
end