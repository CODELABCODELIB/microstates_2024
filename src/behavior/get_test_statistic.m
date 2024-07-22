function [h,p,stats,cutoffs] = get_test_statistic(jid_state_boot,jid_state_orginal, options)
%% Create orginal jid-state
%
% **Usage:**
%   - [h,p,stats] = get_test_statistic(jid_state_boot,jid_state_orginal)
%                 - get_test_statistic(...,'test', 'ttest')
%
% Input(s):
%    jid_state_boot = bootstrapped jid-states 
%           (shape: bin_size*bin_size x number of bootstraps)
%    jid_state_orginal = original jid-state (shape: bin_size*bin_size x 1)
%
% Optional input parameter(s):
%    test (string, default : 'signrank') = chosen test options include
%         'ttest' or 'signrank'
%    alpha (default : 0.05) = alpha for the test
%
% Output(s):
%   h (logical) = result of the hypothesis test
%       h = 1, this indicates rejection of the null hypothesis
%       h = 0, this indicates a failure to reject the null hypothesis
%   p = p-value of the test, returned as a nonnegative scalar from 0 to 1
%   stats (struct) = Test statistics. 
%       In signrank it contains z-value and sign
%       In ttest it contains t-value (tstat), Degrees of freedom (df) and
%       standard deviation (sd)
%
%
% Ruchella Kock, Leiden Unviersity, 12/09/2023
arguments
    jid_state_boot;
    jid_state_orginal
    options.test char = 'prctile';
    options.alpha = 0.05;
    options.thresh (1,2) = [NaN, NaN];
end
if any(isnan([NaN, NaN]))
    options.thresh = [100*(options.alpha/2), 100*(1-(options.alpha/2))];
end
h = NaN; p = NaN; stats = NaN;
% perform ttest
if strcmp(options.test, 'ttest')
    [h,p,~,stats] = ttest(jid_state_boot,jid_state_orginal,'Alpha',options.alpha);
% perform signrank test
elseif strcmp(options.test, 'ttest2')
    [h,p,~,stats] = ttest2(jid_state_boot,jid_state_orginal,'Alpha',options.alpha);
elseif strcmp(options.test, 'signrank')
    [p,h,stats] = signtest(jid_state_boot,jid_state_orginal,'Alpha', options.alpha);
elseif strcmp(options.test, 'prctile')
    cutoffs = prctile(jid_state_boot,options.thresh,2);
    stats = jid_state_orginal<cutoffs(1) || jid_state_orginal>cutoffs(2);
end
end