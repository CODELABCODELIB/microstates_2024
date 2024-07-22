function [jid_state_boot,probabilities_boot] = boot_state_prob(taps,microstates,dt_dt,n_states,options)
arguments
    taps;
    microstates;
    dt_dt;
    n_states;
    options.prob_density logical = 0;
    options.n_boots = 1000;
    options.duration = 1;
    options.pool_method = 'mean';
end
[xi] = create_jid_grid();
probabilities_boot = zeros(n_states, length(taps)-2,options.n_boots);
jid_state_boot = cell(1,n_states);
for boot=1:options.n_boots
    for triad_idx = 1:length(taps)-2
        triad = taps(triad_idx:triad_idx+2);
        duration_triad = (triad(3)-triad(1));
        boot_start = randi(length(microstates)-duration_triad);
        boot_end = boot_start+duration_triad;
        if options.prob_density
            [f,xi] = ksdensity(microstates(boot_start:boot_end), [1:n_states],'Bandwidth', 0.1);
        else
            f = zeros(n_states,1);
            for state=1:n_states
                if (options.duration == 1)
                    f(state) = duration_probability([boot_start, NaN, boot_end],microstates,state);
                elseif options.duration == 2
                    f(state) = occurence_probability([boot_start, NaN, boot_end],microstates,state);
                else
                    f(state) = background_probability([boot_start, NaN, boot_end],microstates,state);
                end
            end
        end
        probabilities_boot(:,triad_idx,boot) = f;
    end
    [dt_dt_tmp,gridx,xi] = assign_tap2bin(dt_dt);
    not_occupied_bins = reshape(~ismember(xi, dt_dt_tmp(:,3:4),'rows'),50,50);
    
    for state=1:n_states
        jid_state_boot_tmp = zeros(50,50);
        [pooled] = pool_duplicates(dt_dt_tmp,probabilities_boot(state,:,boot), 'pool_method',options.pool_method); % pool multiple values in each bi
        for sel_tap=1:size(pooled,1)
            jid_state_boot_tmp(gridx == dt_dt_tmp(sel_tap,3),gridx == dt_dt_tmp(sel_tap,4)) = pooled(sel_tap);
        end
        jid_state_boot_tmp(not_occupied_bins) = NaN;
        jid_state_boot{state}(:,boot) = reshape(jid_state_boot_tmp,2500,1);
    end
end
end