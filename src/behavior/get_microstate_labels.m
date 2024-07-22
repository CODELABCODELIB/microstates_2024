function [microstatestap,dt_dt,n_states,taps] = get_microstate_labels(EEG,indexes,options)
%% Get microstate labels surrounding a tap
%
% **Usage:**
%   - [jid_state] = get_microstate_labels(ALLEEG,pp)
%
% Input(s):
%    EEG = EEG struct
%    indexes cell = latencies where tap sequences take place, if there are
%                   multiple gaps then the cell is of size (1 x number of gaps)
%
% Output(s):
%   microstatestap = labels of microstates around each tap (shape : 1 x number of taps)
%   dt_dt = ITIs K and K+1 (shape : number of taps x 2)
%   n_states = number of microstates 
%
% Ruchella Kock, Leiden Unviersity, 12/09/2023
arguments
    EEG
    indexes
    options.shuffle logical = 0;
end
%% find the taps and ITIs
[taps] = find_taps(EEG, indexes);
[dt_dt,taps] = calculate_ITI_K_ITI_K1(taps, 'shuffle', options.shuffle);

%% select the active microstate at the time of the tap
microstatestap = EEG.microstate.fit.labels(taps);
%% get the number of microstates (regardless if it coincided with a tap)
% all microstates
microstates = EEG.microstate.fit.labels;
% number of states
n_states = max(microstates);
end