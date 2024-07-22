function [Latency_EOStart,Latency_ECStart] = auto_EO_EC_indexer(EEG, labels, latencies, options)
%% Programatically find eyes open or eyes closed indexes based on blink assumptions
%
% **Usage:**
%   - auto_EO_EC_indexer(EEG, labels, latencies)
%   - auto_EO_EC_indexer(..., labels, latencies)
%   - auto_EO_EC_indexer(... , 'Threshold', 5)
%   - auto_EO_EC_indexer(... , 'plot_data', 1)
%
% Input(s):
%    EEG = EEG struct
%
% Optional input parameter(s):
%    labels = Type of events that took place ([EEG.event.type] or [EEG.urevent.type])
%    latencies = Time indexes that all events took place in ms ([EEG.event.latency] or [EEG.urevent.latency])
%    threshold (default 5) = Minimum number of blinks required to accept the index
%            e.g. strict threshold of 0 then there should at some point be absolutely no blinks for 6 minutes.
%            threshold of 5 then there may be a few blinks perhaps participants accidentally opened their eyes a few times
%    plot_data (default 0) = Plot blinks with estimated thresholds
%
% Output(s):
%   Latency_EOStart = Estimated latency start eyes open in ms (if NaN not able to identify)
%   Latency_ECStart = Estimated latency start eyes closed in ms (if NaN not able to identify)
%
% Arko Ghosh, Leiden University, 12/07/2018
% Edited 18/08/2023 Ruchella Kock, Leiden University
arguments
    EEG struct;
    labels cell = [];
    latencies = [];
    options.threshold {mustBeGreaterThanOrEqual(options.threshold,0)} = 5;
    options.plot_data logical = 0;
    options.save_results logical = 1;
    options.save_path = '';
end
%% data checks
% case 2 forgot to send triggers in the right order
% check if there are 64 electrodes, ocular electrodes need to be present
if size(EEG.data,1) < 64
    error('Missing blink electrodes')
end
if isempty(labels)
    labels = {EEG.urevent.type}';
end
if isempty(latencies)
    latencies = [EEG.urevent.latency]';
end
Latency_EOStart = NaN;
Latency_ECStart = NaN;
%% run ICA identify artifacts
try
    EEG.icaquant1 = icablinkmetrics(EEG, 'ArtifactChannel', EEG.data(63,:), 'Alpha', 0.001, 'VisualizeData', 'False');
    EEG.icaquant2 = icablinkmetrics(EEG, 'ArtifactChannel', EEG.data(64,:), 'Alpha', 0.001, 'VisualizeData', 'False');
    Rej_Comp = unique([EEG.icaquant1.identifiedcomponents,EEG.icaquant2.identifiedcomponents]);
    if ~isempty(Rej_Comp)
        EEG = pop_subcomp(EEG,Rej_Comp,0);
    end
catch
    fprintf('Failed getting ICA blink metrics')
    return
end
%% Set Testing_T which truncates the data to only sequences with actual potential resting state data
% resting state always takes place before reaction time session, find resting state start time
labels_RT = strcmp(labels, 'T  1_on');

if sum(labels_RT)>1
    %Tesing_T = latencies(min(find(labels_RT == true))); % Resting state was measured prior to this time
    Testing_T = latencies(find(labels_RT, 1));
    % added!
elseif sum(strcmp(labels, 'S 32')) && sum(strcmp(labels, 'S 36'))
    Testing_T = latencies(find(strcmp(labels, 'S 36'), 1, 'last'));
else
    fprintf('No resting state data')
    return
end
% incase the reaction time trigger was accidentally turned on before the resting state data
% if Testing_T < (6*1000*60) && sum(strcmp(labels, 'S 32')) && sum(strcmp(labels, 'S 36'))
%     Testing_T = latencies(find(strcmp(labels, 'S 36'), 1, 'last'));
% end
%% Identify blinks
Blinks_t = unique([EEG.icaquant1.artifactlatencies,EEG.icaquant2.artifactlatencies]);
% remove blinks that are too close to each other
Blinks_t(diff(Blinks_t)<50)=[];
% Identify blinks before Tesing_T
Blinks_t_under_Testing_T = Blinks_t(Blinks_t<Testing_T);

% find a 6 min window with least amount of blinks and that must be EC
Times_range = 1:Testing_T;
Times_blink = ismember(Times_range, Blinks_t_under_Testing_T);
tt = 0;
% every 1 second check if the next 6 minutes (6*1000*60 = 360000 ms) contains blinks
for t = 1:1000:(Times_range(end)-(6*1000*60))
    period_b = and(Times_range>Times_range(t), Times_range<(Times_range(t)+(6*1000*60)));
    tt = tt+1;
    % number of blinks that takes place over the next 6 minutes
    blink_v(tt) = sum(Times_blink(period_b));
    % time where the blink check took place (every 1 second) saved in ms
    time_v(tt) = Times_range(t);
end
if exist('blink_v','var')
    %% Time eyes closed EC condition started (Select time where almost no blinks took place)
    [smallest_n_blinks, Idx_v] = min(blink_v);
    % added!
    % If smallest number of blinks is below (or equal to) a set threshold then accept the approximation
    if smallest_n_blinks <= options.threshold
        % The presumed latency of eyes closed
        Latency_ECStart = time_v(Idx_v);
    else
        fprintf('Smallest number of blinks(%d) is higher than threshold\n', smallest_n_blinks)
        Latency_ECStart = NaN;
    end
else
    fprintf('Reaction time trigger was accidentally turned on before the resting state data')
    return
end
%% Time eyes open EO condition started
Latency_pre_t_on = latencies<Testing_T;
Latencies_6_min = diff(latencies)>(1000*60*6);

if ~isnan(Latency_ECStart)
    % Check if any of the '0' triggers overlap with the presumed eyes closed period
    % Find resting state events that do not belong to the presumed eyes closed period
    Latency_0 = or((Latency_ECStart>latencies),((Latency_ECStart+(1000*6*60))<latencies));
    % Find resting state events taking place every 6 minutes before reaction time experiment: and(Latency_pre_t_on(1:end-1),Latencies_6_min)
    % that overlap with the 0 triggers: Latency_0(1:end-1)
    Log_latencies_EO = [and(and(Latency_pre_t_on(1:end-1),Latencies_6_min),Latency_0(1:end-1))' false];
    
    % if only one S 32/0 event remains, that is the eyes open start!
    if sum(Log_latencies_EO) == 1
        Latency_EOStart = latencies(Log_latencies_EO);
    elseif sum(Log_latencies_EO) == 2
        % if two events S 32/0 remain, eyes open was probably repeated, so take the last repeat
        Latency_EOStart = max(latencies(Log_latencies_EO));
    else
        fprintf('all efforts to find EO starting latency were unsucessful')
        Latency_EOStart = NaN;
        Latency_ECStart = NaN;
    end
end
%% plot blinks
if options.plot_data
    fig = figure('name', 'Blink pattern');
    scatter(time_v/(1000*60), blink_v)
    xlabel('Time (Minutes)'); ylabel('Number of blinks');
    hold on;
    scatter(Blinks_t_under_Testing_T/(1000*60), ones(1,length(Blinks_t_under_Testing_T)), 'xr');
    yline(options.threshold, '--')
    scatter(latencies(Latencies_6_min)/(1000*60),zeros(length(latencies(Latencies_6_min)),1), '^k', 'filled')
    
    if ~isnan(Latency_EOStart) || ~isnan(Latency_ECStart)
        xline(Latency_EOStart/(1000*60), 'b', 'Linewidth', 5)
        xline(Latency_ECStart/(1000*60), 'r',  'Linewidth', 5)
        legend({'Number of blinks over 6 minutes(calculated every 1 second)','Blinks', 'Threshold','6 minutes interval','EO start', 'EC start'}, 'location', 'southoutside')
    else
        legend({'Number of blinks over 6 minutes(calculated every 1 second)','Blinks', 'Threshold','6 minutes interval'}, 'location', 'southoutside')
    end
    tmp = split(EEG.comments, ':');
    participant = erase(strrep(strcat(tmp{2}(1:6), ';',EEG.filename), '_', '-'), {'.set', ' '});
    title(sprintf('Blinks during resting state sessions EO and EC: Sub - %s',participant))
    if (options.save_results)
        saveas(fig, strcat(options.save_path,sprintf('blink_resting_state_%s.svg',participant)))
    end
end