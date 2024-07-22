function [EEG,resting_dat_present] = EO_EC_recomender(EEG)
%% Select the start times for resting state eyes open or closed either via 
%   either the manual trigger input or automatic blink analysis 
%
% **Usage:**
%   - EO_EC_recomender(EEG)
%
% Input(s):
%    EEG = EEG struct
%
% Output(s):
%    EEG = EEG struct with resting state conditions start times
%    resting_dat_present = logical is there any resting state data
%
% Requires:
%    manual_EO_EC_indexer.m
%    auto_EO_EC_indexer.m
%
% Quick cheatsheet manual triggers:
% Type: S 33/1 EO
% Type: S 34/2 EC
% Type: S 32/0 start EO or EC
% Type: S 36/4 end EO or EC
% possible example sequences
%     - EO 1-0-4 EC 2-0-4
%     - EC 2-0-4 EO 1-0-4
%     - EO 1-0-4 1-0-4  EC 2-0-4
%
% Arko Ghosh, Leiden University, 12/07/2018
% Edited 18/08/2023 Ruchella Kock, Leiden University
%%
labels = {EEG.event.type}';
latencies = [EEG.event.latency]';
labels_1 = strcmp(labels, 'S 33');
labels_2 = strcmp(labels, 'S 34');
labels_0 = strcmp(labels, 'S 32');
labels_4 = strcmp(labels, 'S 36');

% first check if there is any resting state data in the file if not return
resting_dat_present = 1;
if ~sum(labels_0) 
    resting_dat_present = 0;
    fprintf('No resting state data present')
    return
end

[EEG.Resting_Trigger_EOStart,EEG.Resting_Trigger_ECStart] = manual_EO_EC_indexer(latencies, labels_0,labels_1,labels_2,labels_4);
[EEG.Resting_Blink_EOStart,EEG.Resting_Blink_ECStart] = auto_EO_EC_indexer(EEG, labels, latencies, 'plot_data', 1, 'save_path', '/home/ruchella/microstates_2023/results/figures/');

% Recomend Resting stated decisions
% most accurate is the blink based start times
if ~isnan(EEG.Resting_Blink_EOStart)
    EEG.Reco_EO_start = EEG.Resting_Blink_EOStart;
% otherwise select the manual trigger start times
elseif ~isnan(EEG.Resting_Trigger_EOStart)
    EEG.Reco_EO_start = EEG.Resting_Trigger_EOStart;
else
    EEG.Reco_EO_start = NaN;
end

if ~isnan(EEG.Resting_Blink_ECStart)
    EEG.Reco_EC_start = EEG.Resting_Blink_ECStart;
elseif ~isnan(EEG.Resting_Trigger_ECStart)
    EEG.Reco_EC_start = EEG.Resting_Trigger_ECStart;
else
    EEG.Reco_EC_start = NaN;
end
end