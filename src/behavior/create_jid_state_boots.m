function [jid_boot] = create_jid_state_boots(EEG,indexes,state,options)
%% Create bootstrapped jid-states
%
% **Usage:**
%   - [jid_boot] = create_boots(labels,dt_dt,state)
%                - create_boots(...,'shuffle', 0)
%
% Input(s):
%    labels = microstate labels that get randomly sampled
%    dt_dt = ITIs K and K+1 (shape : number of taps x 2)
%    state = double given microstate
%
% Optional input parameter(s):
%   n_boots = number of bootstraps to create
%   bin_size = bin size of the jid-state
%   reshape (default : 0) = 1 reshape each boot to size bin_size x bin_size, 0 otherwise
%
% Output(s):
%   jid_boot = bootstrapped jid-states 
%       (shape: bin_size*bin_size x n_boots) - if reshaped is 0 
%       (shape: bin_size x bin_size x n_boots) - if reshaped is 1
%
% Requires:
%   create_jid_grid.m
%
% Ruchella Kock, Leiden Unviersity, 12/09/2023
arguments
    EEG;
    indexes;
    state;
    options.n_boots = 1000;
    options.bin_size = 50;
    options.reshape logical = 0;
end
%% initialize data
jid_boot = zeros(options.bin_size*options.bin_size,options.n_boots);
[xi] = create_jid_grid(options.bin_size);
%% create bootstrapped jid-states
for boot=1:options.n_boots
    % get indexes of the randomly selected triads
    [labels,dt_dt,n_states] = get_microstate_labels(EEG,indexes, 'shuffle', 1);
    % get the ITI's K and K+1 that belong to each state
    selected_triads = dt_dt(labels == state,:);
    % there needs to be atleast 1 per microstate to make a probability plot
    % if its less than 1 then remove that microstate
    if size(selected_triads,1) <= 1
        warning('JID-state %d = 0\n',state)
        jid_boot(:)= deal(NaN);
        return;
    else
        % per microstate create the JID-state with the bootstrapped data
        jid_boot(:,boot) = ksdensity(selected_triads,xi, 'Bandwidth', 0.1);
    end
    if options.reshape
        jid_boot(:,boot) = reshape(jid_boot(:,boot),options.bin_size,options.bin_size);
    end    
end