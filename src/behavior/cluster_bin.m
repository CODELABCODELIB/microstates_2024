function [stable_basis, stable_loading,reconstruct] = cluster_bin(ALLEEG,labels,pps_files_in_cluster, options)
% Input(s):
%    ALLEEG struct = all EEG structs for all participants with microstate field
%           (shape : number of participants x 2) where the first column
%           contains the participant number (unused in function) and second
%           the EEG struct for each participant WITH microstates field
%
% Optional input parameter(s):
%   save_results logical (default : 1) = 1 save each jid-state and boot
%           jid-states, save the h,p and stats, 0 otherwise
%   save_path (default : current folder) = path to save results
%   start_run (default : 1) = index to start with in ALLEEG
%   end_run (default : size(ALLEEG,1)) = index to end with in ALLEEG
%   repetitions (default : 100) = number of times to repeat the nnmf pipeline
%   plot_res logical (default : 1) = 1 plot nnmf results, 0 otherwise
%   parameter (default : 'tp') = parameter to calculate 'duration' transition duration
%            or 'tp' transition probability
%
% % Output(s):
%
% Requires
%   assign_tap2bin.m
%   calculate_TP.m
%   nnmf_pipeline_spams
%
% Author(s)
% Ruchella Kock, Leiden University, 15/11/2023

arguments
    ALLEEG;
    labels = [];
    pps_files_in_cluster = [];
    options.save_results logical = 1;
    options.save_path char = '.';
    options.file = '';
    options.repetitions = 50;
    options.start_run = 1;
    options.end_run = size(ALLEEG,1);
    options.plot_res logical = 0;
    options.parameter = 'tp';
    options.threshold = 0;
    options.n_bins = 50;
    options.remove_states = 0;
    options.noisy_states = [2,3,9,4,7];
end

for pp=options.start_run:options.end_run
%     stable_basis = [];   stable_loading = [];  reconstruct = [];
    % select one participant EEG
    EEG = ALLEEG{pp,2};
    % find tap indexes
    [taps] = find_taps(EEG, ALLEEG{pp,3});
    % identify the jid bins
    [dt_dt,~] = calculate_ITI_K_ITI_K1(taps, 'shuffle', 0);
    % assign taps to the jid bins
    [dt_dt,gridx,xi] = assign_tap2bin(dt_dt,'BINS',options.n_bins);
    % identify bins with no behavior
    %     not_occupied_bins = reshape(~ismember(xi, dt_dt(:,3:4),'rows'),50,50);
    %% populate bins with microstate sequences
    jid_microstates = cell(options.n_bins,options.n_bins);
    % get microstate labels
    microstates = EEG.microstate.fit.labels;
    
    % reorder microstates based on pop clusters
    if ~isempty(labels) && ~isempty(pps_files_in_cluster)
        EEG.microstate.pop_labels = labels(contains(pps_files_in_cluster,EEG.filepath(end-3:end)));
        [indices,reordered_microstates_labels] = reorder_microstates(EEG.microstate.pop_labels);
        [microstates] = replace_microstate_labels(microstates,indices);
        EEG.microstate.reordered_indices = indices;
        EEG.microstate.reordered_microstates_labels = reordered_microstates_labels;
        EEG.microstate.reordered_states = microstates;
        EEG.microstate.filename = EEG.filepath(end-3:end);
    end
    
    % populate jid_state_original_tmp with the microstate sequences during each triad
    % -2 for the last 2 taps that do not fully form a triad
    for triad_idx = 1:length(taps)-2
        triad = taps(triad_idx:triad_idx+2); % triad indexes
        % select microstate sequences during the triad
        mstate_sequence = microstates(triad(1):triad(3));
        % check if the bin is occupied already if not add the sequences
        if isempty(jid_microstates{gridx == dt_dt(triad_idx,3),gridx == dt_dt(triad_idx,4)})
            jid_microstates{gridx == dt_dt(triad_idx,3),gridx == dt_dt(triad_idx,4)} = {mstate_sequence};
            % if it is occupied concat the new sequence to existing one(s)
        else
            jid_microstates{gridx == dt_dt(triad_idx,3),gridx == dt_dt(triad_idx,4)} = cat(1,jid_microstates{gridx == dt_dt(triad_idx,3),gridx == dt_dt(triad_idx,4)},{mstate_sequence});
        end
    end
    % flatten
    jid_microstates = reshape(jid_microstates, options.n_bins*options.n_bins,1);
    %% get transition probabilities for each bin
    % nan of number bins x total number of transition probabilities
    % there are N_states prototypes so the combinates are N_states*N_states
    N_states_original = size(EEG.microstate.prototypes,2);
    if options.remove_states && ~isempty(labels)
        states_2_remove = ismember(reordered_microstates_labels,options.noisy_states);
        EEG.microstate.removed_states = states_2_remove;
        N_states = sum(~states_2_remove);
    else
        N_states = N_states_original;
    end
    TP_jid = nan(options.n_bins*options.n_bins,N_states*N_states);
    TP_tmp_jid = nan(options.n_bins*options.n_bins,N_states*N_states);
    for bin=1:options.n_bins*options.n_bins
        binseq = jid_microstates{bin};
        if ~isempty(binseq)
            % concat all sequences in one bin following each other and
            % calculate the transition probabilities
            [TP] = calculate_TP(N_states_original,[binseq{:}],'parameter',options.parameter);
            if options.remove_states && ~isempty(labels)
                TP = TP(~states_2_remove,~states_2_remove);
            end
            TP_jid(bin,:) = reshape(TP, N_states*N_states,1);
            
            % threshold based on empty values easier to check with duration
            % because these are nans wheras TP are set to zero if bin is empty
            if strcmp(options.parameter, 'tp')
                [TP_tmp] = calculate_TP(N_states_original,[binseq{:}],'parameter','duration');
                if options.remove_states && ~isempty(labels)
                    TP_tmp = TP_tmp(~states_2_remove,~states_2_remove);
                end
                TP_tmp_jid(bin,:) = reshape(TP_tmp, N_states*N_states,1);
            end
        end
    end
    if contains(options.parameter, 'duration')
        TP_tmp_jid = TP_jid;
    end
    if options.save_results
        save(sprintf('%s/%s_jid_microstates_%d',options.save_path,options.file, pp),'jid_microstates')
        save(sprintf('%s/%s_TP_jid_%d',options.save_path,options.file,pp),'TP_jid')
    end
    % remove the completely empty bins before NNMF
    removed_bins = isnan(TP_tmp_jid);
    kept_bins_idx = find(~all(removed_bins'));
    TP_jid = TP_jid(kept_bins_idx,:);
    TP_tmp_jid  = TP_tmp_jid(kept_bins_idx,:);
    
    if options.threshold
        % remove the bins that are more than 75% empty compared to the other bins
        num_empty_cells = sum(isnan(TP_tmp_jid),2);
        kept_bins_2 = num_empty_cells < quantile(num_empty_cells, options.threshold);
        TP_jid = TP_jid(kept_bins_2,:);
        
        % index of the kept bins to assign later
        kept_bins = logical(zeros(options.n_bins*options.n_bins,1));
        kept_bins(kept_bins_idx(kept_bins_2)) = 1;
    else
        kept_bins = logical(zeros(options.n_bins*options.n_bins,1));
        kept_bins(kept_bins_idx) = 1;
    end
    
    if contains(options.parameter, 'duration') && ~isempty(TP_jid)
        TP_jid(isnan(TP_jid)) = 0.000000000000000001;
        if min(log10(TP_jid), [],'all') < 0
            TP_jid = log10(TP_jid) + abs(min(log10(TP_jid), [],'all'));
        else
            TP_jid = log10(TP_jid);
        end
    end
    %% NNMF
    basis_all = {};
    loadings_all = {};
    stable_basis = [];
    stable_loading = [];
    reconstruct = [];
    if ~isempty(TP_jid)
        [~, ~, ~, test_err] = nnmf_cv(TP_jid, 'repetitions', options.repetitions);
        [best_k_overall]  = choose_best_k({test_err}, 1);
        % perform nnmf multiple times
        for rep = 1:100
            [W, H] = perform_nnmf(TP_jid, best_k_overall);
            basis_all{rep,1} = W;
            loadings_all{rep,1} = H;
        end
        [stable_basis, stable_loading] = stable_nnmf(basis_all,loadings_all, 1);
        stable_loading = full(stable_loading);
        %% reconstruct the jid bins
        reconstruct = nan(options.n_bins*options.n_bins,best_k_overall);
        for k=1:best_k_overall
            reconstruct(kept_bins,k) = stable_basis(:,k);
        end
    end
    if options.save_results
        save(sprintf('%s/%s_stable_basis_%d',options.save_path,options.file, pp),'stable_basis')
        save(sprintf('%s/%s_stable_loading_%d',options.save_path,options.file,pp),'stable_loading')
        save(sprintf('%s/%s_reconstruct_%d',options.save_path,options.file,pp),'reconstruct')
        eeg_microstate = EEG.microstate;
        save(sprintf('%s/%s_microstate_%d',options.save_path,options.file,pp),'eeg_microstate')
    end
    %% display NNMF results
    if options.plot_res && ~isempty(stable_loading)
        figure;
        tiledlayout(2,best_k_overall);
        for k=1:best_k_overall
            nexttile;
            imagesc(reshape(stable_loading(k,:),N_states,N_states)); colorbar;
            %         caxis([0,0.4])
            set(gca, 'YDir', 'normal')
            title(sprintf('rank %d meta-transitions',k))
        end
        for k=1:best_k_overall
            nexttile;
            jid = reshape(reconstruct(:,k),options.n_bins,options.n_bins);
            plot_jid(jid); colorbar;
            title(sprintf('rank %d meta-jidtransitions',k))
            %         caxis([0,1.5])
            sgtitle(sprintf('Sub %d',pp));
        end
    end
end