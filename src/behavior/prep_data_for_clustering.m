function [parameters,jid_state_params] = prep_data_for_clustering(ALLEEG,stats,boots,original,options)
%% Prepare data for clustering
%
% **Usage:**
%   - prep_data_for_clustering(ALLEEG)
%   - prep_data_for_clustering(...,stats)
%
% Input(s):
%    ALLEEG = EEG structs of all participants with microstate field
%
% Optional input parameter(s):
%   stats cell = statistics results, such as FDR corrected jid-states - shape(number of pps x 1)
%   boots cell = jid-states bootstraps - shape(number of pps x 1)
%   original cell = jid-state probabilities - shape(number of pps x 1)
%
% Output(s):
%   parameters.all_prototypes = all the prototypes;
%   parameters.pps_in_cluster = which participant the prototypes belong to;
%   parameters.all_jids = jid-states;
%   parameters.all_selected = ITIs used to create the JID-state;
%   parameters.Duration = Duration microstate parameter;
%   parameters.Occurence = Occurence microstate parameter;
%   parameters.Coverage = Coverage microstate parameter;
%   parameters.GEV = GEV  microstate parameter;
%   parameters.all_stats = statistics results;
%   parameters.all_boots = jid-states bootstraps;
%   parameters.all_originals = jid-state probabilities;
%
% Author: R.M.D. Kock, Leiden University, 04/12/2023

arguments
    ALLEEG;
    stats cell = {};
    boots cell = {};
    original cell= {};
    options.parameter = [];
    options.save_results logical = 0;
    options.save_path ='.';
    options.file = '';
    options.catdim = 3;
end
%% prepare the variables
all_prototypes = []; pps_in_cluster = [];
all_jids = []; all_selected = [];
all_stats = []; all_boots = []; all_originals = [];
Duration = []; Occurence = []; Coverage = []; GEV = [];
count = 1;
[jid_state_params] = create_jid_states(ALLEEG);
for pp =1:size(ALLEEG,1)
    EEG = ALLEEG{pp,2};
    n_states = max(EEG.microstate.fit.labels);   
    has_jid = jid_state_params{pp,5};
    
    %% JID state
    % all_jids contains the jid states
    all_jids = cat(options.catdim,all_jids,jid_state_params{pp,3}(has_jid));
    % all_selected contains the ITIs used to create the JID-state
    all_selected = cat(options.catdim,all_selected,jid_state_params{pp,4}(has_jid));
    %% prototypes from each participant to cluster at the population level
    % select the prototype
    tmp_cat = EEG.microstate.Res.A_all{n_states-1, 1}; 
    % remove the prototypes of the microstates that are absent
    tmp_cat = tmp_cat(:, has_jid);
    % all_clusters gives all the available prototypes
    all_prototypes = cat(options.catdim,all_prototypes,tmp_cat);
    % pps_in_cluster keeps track of which participant the prototypes
    % belonged to
    pps_in_cluster = cat(options.catdim,pps_in_cluster,ones(1,size(tmp_cat,2))*pp);
    %% microstate parameters
    Duration = cat(options.catdim,Duration,EEG.microstate.stats.Duration(has_jid));
    Occurence = cat(options.catdim,Occurence,EEG.microstate.stats.Occurence(has_jid));
    Coverage = cat(options.catdim,Coverage,EEG.microstate.stats.Coverage(has_jid));
    GEV = cat(options.catdim,GEV,EEG.microstate.stats.GEV(has_jid));
    %% extra parameters
    if ~isempty(stats)
        all_stats = cat(options.catdim,all_stats,stats{pp}(has_jid));
    end
    if ~isempty(boots)
        all_boots = cat(options.catdim,all_boots,boots{pp}(has_jid));
    end
    if ~isempty(original)
        all_originals = cat(options.catdim,all_originals,original{pp}(has_jid));
    end
    count = count+1;
end
parameters = struct();

parameters.all_prototypes = all_prototypes;
parameters.pps_in_cluster = pps_in_cluster;
parameters.all_jids = all_jids;
parameters.all_selected = all_selected;
parameters.Duration = Duration;
parameters.Occurence = Occurence;
parameters.Coverage = Coverage;
parameters.GEV = GEV;
parameters.all_stats = all_stats;
parameters.all_boots = all_boots;
parameters.all_originals = all_originals;
parameters.jid_state_params = jid_state_params;
%% save results
if options.save_results
    save(sprintf('%s/%s',options.save_path,options.file), 'parameters')
end
end