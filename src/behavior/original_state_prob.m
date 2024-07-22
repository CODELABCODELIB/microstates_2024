function [jid_state_original,probabilities] = original_state_prob(taps,microstates,dt_dt,n_states,options)
arguments
    taps;
    microstates;
    dt_dt;
    n_states;
    options.prob_density logical = 0;
    options.duration = 1;
    options.pool_method = 'mean';
end
[xi] = create_jid_grid();
probabilities = zeros(n_states, length(taps)-2);
for triad_idx = 1:length(taps)-2
    triad = taps(triad_idx:triad_idx+2);
    if options.prob_density
        [f] = ksdensity(microstates(triad(1):triad(3)), [1:n_states],'Bandwidth', 0.1);
    else
        f = zeros(n_states,1);
        for state=1:n_states
            if (options.duration == 1)
                f(state) = duration_probability(triad,microstates,state);
                disp('duration')
            elseif (options.duration == 2)
                f(state) = occurence_probability(triad,microstates,state);
                disp('occurence')
            elseif (options.duration == 3)
                f(state) = coverage_probability(triad,microstates,state);
                disp('coverage')
            else
                f(state) = background_probability(triad,microstates,state);
            end
        end
    end
    probabilities(:,triad_idx) = f;
end
%% pooled probabilities
[dt_dt,gridx,xi] = assign_tap2bin(dt_dt);
not_occupied_bins = reshape(~ismember(xi, dt_dt(:,3:4),'rows'),50,50);
jid_state_original = cell(1,n_states);
for state=1:n_states
    jid_state_original_tmp = zeros(50,50);
    [pooled] = pool_duplicates(dt_dt,probabilities(state,:), 'pool_method',options.pool_method); % pool multiple values in each bi
    for sel_tap=1:size(pooled,1)
        jid_state_original_tmp(gridx == dt_dt(sel_tap,3),gridx == dt_dt(sel_tap,4)) = pooled(sel_tap);
    end
    jid_state_original_tmp(not_occupied_bins) = NaN;
    jid_state_original{state} = reshape(jid_state_original_tmp,2500,1);
end
end