function [Resting_Trigger_EOStart,Resting_Trigger_ECStart] = manual_EO_EC_indexer(latencies, labels_0,labels_1,labels_2,labels_4)
%% Identify resting state manually added start events 
%
% **Usage:**
%   - trigger_EO_EC_indexes(latencies, labels_0,labels_1,labels_2,labels_4)
%
% Input(s):
%   latencies = Time indexes that all events took place in ms ([EEG.event.latency] or [EEG.urevent.latency])
%   labels_0 = Logical array of event Type: S 32/0 start EO or EC
%   labels_1 = Logical array of event Type: S 33/1 EO taking place
%   labels_2 = Logical array of event Type: S S 34/2 EC taking place
%   labels_4 = Logical array of event Type: S 36/4 end EO or EC
%
% Output(s):
%   Resting_Trigger_EOStart = Reported latency start eyes open in ms (if NaN not able to identify)
%   Resting_Trigger_ECStart = Reported latency start eyes closed in ms (if NaN not able to identify)
%
% Quick cheatsheet:
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
%% Check for EO - 1-0-4 sequences
possible_open = and(and(labels_1(1:end-2),labels_0(2:end-1)),labels_4(3:end));
% find the index of the EO condition
Idx = find(possible_open == true);
% check if there isnt a stray trigger between the indices
if isempty(Idx)
    possible_open = and(and(labels_1(1:end-3),labels_0(2:end-2)),labels_4(4:end));
    Idx = find(possible_open == true);
end
% possible duration of the session is 6 minutes
% Calculate and convert dur in minutes (latencies in ms / 1000 * 60)
% Idx + 2 select the latency of the end of the session
% Idx + 1 select the latency of the start of the session
dur = (latencies(Idx+2)-latencies(Idx+1))./(1000*60);
% if there is only one start then select that as the resting trigger EO start
% if there are more than one then the last repeated session that lasted 6 minutes is selected as the start
if and(dur >= 6, length(Idx) == 1)
    Resting_Trigger_EOStart = latencies(Idx+1);
    % added!
elseif and(dur >= 6, length(Idx) > 1)
    Resting_Trigger_EOStart = latencies(Idx(end)+1);
else
    Resting_Trigger_EOStart = NaN;
end
%% Check for EC - 2-0-4 sequences
possible_closed = and(and(labels_2(1:end-2),labels_0(2:end-1)),labels_4(3:end));
Idx_1 = find(possible_closed == true);

dur = (latencies(Idx_1+2)-latencies(Idx_1+1))./(1000*60);
if and(dur >= 6, length(Idx_1) == 1)
    Resting_Trigger_ECStart = latencies(Idx_1+1);
    % added!
elseif and(dur >= 6, length(Idx) > 1)
    Resting_Trigger_ECStart = latencies(Idx(end)+1);
else
    Resting_Trigger_ECStart = NaN;
end
end