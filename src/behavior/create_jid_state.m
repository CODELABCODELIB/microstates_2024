function [jid_state] = create_jid_state(labels,dt_dt,state,options)
%% Create orginal jid-state
%
% **Usage:**
%   - [jid_state] = create_jid_state(labels,dt_dt,state)
%                 - create_jid_state(...,'bin_size', 40)
%
% Input(s):
%    labels = microstate labels that get randomly sampled
%    dt_dt = ITIs K and K+1 (shape : number of taps x 2)
%    state = double given microstate
%
% Optional input parameter(s):
%    bin_size (default : 50) = bin size of the jid-state 
%    reshape (default : 0) = 1 reshape the data to size bin_size x bin_size, 0 otherwise
%
% Output(s):
%   jid_state = original jid-state 
%       (shape: bin_size*bin_size x 1) - if reshaped is 0 
%       (shape: bin_size x bin_size) - if reshaped is 1
%
% Requires:
%   create_jid_grid.m
%
% Ruchella Kock, Leiden Unviersity, 12/09/2023
arguments
    labels;
    dt_dt;
    state;
    options.bin_size = 50;
    options.reshape logical = 0;
end
[xi] = create_jid_grid(options.bin_size);
% get the ITI's K and K+1 that belong to each state
selected_triads = dt_dt(labels == state,:);
% there needs to be atleast 1 per microstate to make a probability plot
% if its less than 1 then remove that microstate
if size(selected_triads,1) <= 1
    warning('JID-state %d = 0\n',state)
    jid_state(:)= NaN(2500,1);
else
    % per microstate create the JID-state grid
    jid_state = ksdensity(selected_triads,xi, 'Bandwidth', 0.1);
end
% reshape the input to size options.bin_size x options.bin_size
if options.reshape
    jid_state = reshape(jid_state,options.bin_size,options.bin_size);
end
end