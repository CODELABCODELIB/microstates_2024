function [stats_adj, stats_og] = FDR(jid_state_boot, jid_state_original, state,options)
arguments
    jid_state_boot;
    jid_state_original;
    state;
    options.plot_res = 0;
end
%% load the boots and original data

grid_prc = nan(2500,1);
pvals = nan(2500,1);
for bin=1:2500
    % calculate the percentile for the 'real' value based on the distribution
    grid_prc(bin) = invprctile(jid_state_boot{state}(bin,:), jid_state_original{state}(bin));
    % convert to 'p vals'
    if grid_prc(bin)>= 50
        pvals(bin) = ((100-grid_prc(bin))/100)*2;
    else
        pvals(bin) = (grid_prc(bin)/100)*2;
    end
end
% perform fdr on whole grid
fdr_res = mafdr(pvals, 'BHFDR', 1);
% FDR corrected percentage
adj_grid_prc = nan(2500,1);
adj_grid = nan(2500,1);
stats_adj = nan(2500,1);
% remove?
stats_og = nan(2500,1);
for bin=1:2500
    % convert the fdr adjusted pvals to percentages
    if grid_prc(bin)>= 50
        adj_grid_prc(bin) = (1 - fdr_res(bin) / 2) * 100;
    else
        adj_grid_prc(bin) = fdr_res(bin) / 2 * 100;
    end
    % now from the adjusted percentage get back to the real value in the distribution
    adj_grid(bin) = prctile(jid_state_boot{state}(bin,:),adj_grid_prc(bin),2);
    [~,~,stats_adj(bin), cutoffs] = get_test_statistic(jid_state_boot{state}(bin,:),adj_grid(bin),'test', 'prctile');
    [~,~,stats_og(bin), ~] = get_test_statistic(jid_state_boot{state}(bin,:),jid_state_original{state}(bin),'test', 'prctile');
    
    if options.plot_res && stats_adj(bin)
        figure; histogram(jid_state_boot{state}(bin,:))
        xline(jid_state_original{state}(bin), 'LineWidth', 2)
        xline(adj_grid(bin), 'b', 'LineWidth', 2)
        xline(cutoffs(1), 'r')
        xline(cutoffs(2), 'r')
        title(sprintf('bin %d',bin))
    end
end
if options.plot_res
    [xi] = create_jid_grid();
    figure;
    tiledlayout(1,2)
    nexttile;
    plot_jid(xi,reshape(stats_adj,50,50))
    nexttile;
    plot_jid(xi,reshape(stats_og,50,50))
end
