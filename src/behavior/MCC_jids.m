function [ClusteredOut] = MCC_jids(jid_state_t,jid_state_boot, options)
%% Perform spatial multiple comparison correction on the jid state statistics
%
% **Usage:**
%   - [ClusteredOut] = MCC_jids(jid_state_t,jid_state_boot)
%          - MCC_jids(...,'thresh', 95)
%
% Input(s):
%    jid_state_t = real dataset of test statistics (size 2500 x 1)
%    jid_state_boot = bootstrapped dataset (size 2500 x number of bootstraps)
%
% Optional input parameter(s):
%    pval (default : 0.05) = p value threshold for clustering
%    minchan (default : 5) = minimum number of neighbouring pixels needed for clustering
%    thresh (default : 97.5) = bootsrap threshold (considering two-tailed stats)
%    n_boots (default : calculated based of jid_state_boot)= number of bootstraps to consider
%    plot_res (default : 0) = create plot of the real and selected clusters
%
% Output(s):
%   ClusteredOut = selected clusters based on the multiple comparison correction
%
% Requirements:
%   - find_adj_matrix https://github.com/CODELABLEIDEN/TapDataAnalysis
%   - limo_findcluster https://github.com/LIMO-EEG-Toolbox/limo_tools
%
% Author Arko Ghosh, Leiden University
% Edited  12/09/2023 for jid-states, Ruchella Kock, Leiden University
%
arguments
    jid_state_t;
    jid_state_boot;
    options.pval = 0.05;
    options.minchan = 5;
    options.thresh = [];
    options.n_boots = size(jid_state_boot,2);
    options.plot_res logical = 1;
end
% Set defaults
if isempty(options.thresh)
    options.thresh = 100*(1-(options.pval/2));
end
nM = find_adj_matrix(50,1); % Neighboorhood matrix for LIMO
if options.plot_res
    [xi] = create_jid_grid();
end
%% Step 1 multiple comparison correction
% jid_state_boot = jid_boots{pp}';
% jid_state_t = squeeze(t_values(:,state,pp));
% booleanval = jid_state_t>prctile(jid_state_boot,options.thresh,2);
% run limo clustering (2D) for real data
[clusteredVal, numCluster] = limo_findcluster(jid_state_t,nM,options.minchan);
%% Step 2 multiple comparison correction
for boot = 1:options.n_boots
    % compare each boot with the rest of the boots
    tmp_bootrealval = jid_state_boot(:,boot);
%     tmp_bootbootval = jid_state_boot(:,[1:options.n_boots]~=boot);
%     tmp_booleanval = tmp_bootrealval>prctile(tmp_bootbootval,options.thresh,2);
    % run limo clustering (2D) for bootstrapped data
    [tmp_clusterslabel, nclust] = limo_findcluster(tmp_bootrealval,nM,options.minchan);
    
    % get the cluster size of each cluster
    tmp_Dval = [];
    for n = 1:nclust
        tmp_idx = [tmp_clusterslabel==n]; % select the cluster
        tmp_Dval = [tmp_Dval sum(tmp_idx)]; % get the cluster size
    end
    % Collect the boot values
    if ~isempty(tmp_Dval)
        D_boot(boot) = max(tmp_Dval);  % select the largest cluster size
    else
        D_boot(boot) = 0;
    end
end
%% now keep only those clusters which are both real and < pval on the boot based on max D values
if numCluster
    for n = 1:numCluster
        %     idx = [];
        idx = [clusteredVal==n];
        Dval = [sum(jid_state_t(idx))];
        num_chance(n,1) = Dval > prctile(D_boot,options.thresh);
    end
    
    fcluster = [1:numCluster];
    fcluster(~num_chance) = []; % list clusters that survived
    clusteredVal_select = ismember(clusteredVal,fcluster);
    ClusteredOut = reshape(clusteredVal_select,50,50,size(clusteredVal_select,2));
else
    ClusteredOut = NaN(50,50);
end
%%
if options.plot_res
    figure; tiledlayout(1,3);
    nexttile; plot_jid(xi, reshape(jid_state_t,50,50)); colorbar;
    nexttile; plot_jid(xi, reshape(clusteredVal,50,50)); colorbar;
    nexttile; plot_jid(xi, ClusteredOut); colorbar;
end