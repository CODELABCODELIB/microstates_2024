function [stable_basis, stable_loading,reconstruct,TP_jid] = parameters_nnmf(ALLEEG,labels,pps_files_in_cluster, options)
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
    options.threshold = 1;
    options.n_bins = 50;
    options.remove_states = 0;
    options.noisy_states = [2,3,9,4,7];
end

for pp=options.start_run:options.end_run
    stable_basis = [];   stable_loading = [];  reconstruct = [];
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
    
    if ~isempty(labels) && ~isempty(pps_files_in_cluster)
        EEG.microstate.pop_labels = labels(contains(pps_files_in_cluster,EEG.filepath(end-3:end)));
        [indices,reordered_microstates_labels] = reorder_microstates(EEG.microstate.pop_labels);
        [microstates] = replace_microstate_labels(microstates,indices);
        EEG.microstate.reordered_indices = indices;
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
    N_params = 3;
    
    if options.remove_states && ~isempty(labels)
        % find any microstate that was clustered in the population as noisy
        states_2_remove = ismember(reordered_microstates_labels,options.noisy_states);
        EEG.microstate.removed_states = states_2_remove;
        N_states = sum(~states_2_remove);
    else
        N_states = N_states_original;
    end
    
    TP_jid = nan(options.n_bins*options.n_bins,N_states*N_params);
    for bin=1:options.n_bins*options.n_bins
        binseq = jid_microstates{bin};
        if ~isempty(binseq)
            m_param = nan(N_states_original,N_params);
            for state=1:N_states_original
                [runvalue, runs] = my_RLE([binseq{:}]);
                %% 
                if ~isnan(nanmean(runs(runvalue == state)))
                    % Mean Duration
                    m_param(state,1) =  nanmean(runs(runvalue == state)) .* (1000 / EEG.srate);
                    % Occurence
                    m_param(state,2) =  length(runs(runvalue == state))./length([binseq{:}]).* EEG.srate;
                    % Coverage
                    m_param(state,3) = (m_param(state,1) .* m_param(state,2))./ 1000;
                end
            end
            % remove noisy microstates
            if options.remove_states && ~isempty(labels)
                m_param = m_param(~states_2_remove,:);
            end
            
            TP_jid(bin,:) = reshape(m_param, N_states*N_params,1);
        end
    end
    if options.save_results
        save(sprintf('%s/%s_jid_microstates_%d',options.save_path,options.file, pp),'jid_microstates')
        save(sprintf('%s/%s_TP_jid_%d',options.save_path,options.file,pp),'TP_jid')
    end
    %% remove the completely empty bins before NNMF
    removed_bins = isnan(TP_jid);
    kept_bins_idx = find(~all(removed_bins'));
    TP_jid = TP_jid(kept_bins_idx,:);
    
    if options.threshold
        % remove the bins that are more than 75% empty compared to the other bins
        num_empty_cells = sum(isnan(TP_jid),2);
        kept_bins_2 = num_empty_cells < options.threshold ;
        TP_jid = TP_jid(kept_bins_2,:);
        
        % index of the kept bins to assign later
        kept_bins = logical(zeros(options.n_bins*options.n_bins,1));
        kept_bins(kept_bins_idx(kept_bins_2)) = 1;
    else
        kept_bins = logical(zeros(options.n_bins*options.n_bins,1));
        kept_bins(kept_bins_idx) = 1;
    end
    
    
    %     if contains(options.parameter, 'duration')
    %         TP_jid(isnan(TP_jid)) = 0.000000000000000001;
    %         if min(log10(TP_jid), [],'all') < 0
    %             TP_jid = log10(TP_jid) + abs(min(log10(TP_jid), [],'all'));
    %         else
    %             TP_jid = log10(TP_jid);
    %         end
    %     end
    %% NNMF
    basis_all = {};
    loadings_all = {};
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
    %%
    if options.plot_res
        figure;
        tiledlayout(best_k_overall,N_params+1);
        
        for k=1:best_k_overall
            reshaped_stable_loading = reshape(stable_loading(k,:),N_states,N_params);
            nexttile;
            plot(reshaped_stable_loading(:,1)');
            set(gca, 'YDir', 'normal')
            title(sprintf('rank %d meta-durations',k))
            axis square
            
            nexttile;
            plot(reshaped_stable_loading(:,2)');
            set(gca, 'YDir', 'normal')
            title(sprintf('rank %d meta-occurence',k))
            axis square
            
            %             nexttile;
            %             plot(reshaped_stable_loading(:,3)');
            %             set(gca, 'YDir', 'normal')
            %             title(sprintf('rank %d meta-coverage',k))
            %             axis square
            
            nexttile;
            jid = reshape(reconstruct(:,k),options.n_bins,options.n_bins);
            plot_jid(jid, 'BINS', options.n_bins); colorbar;
            colorbar;
            title(sprintf('rank %d meta-jidtransitions',k))
            caxis([0,1])
            sgtitle(sprintf('Sub %d',pp));
            axis square
        end
    end
end