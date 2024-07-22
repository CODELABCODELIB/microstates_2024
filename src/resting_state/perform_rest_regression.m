function [results,params_per_pps] = perform_rest_regression(data_path,type,AllJID,rest_subjects,options)
%% Perform regression between JID 1 month before measurement and microstate parameter
%
% **Usage:** [params_per_pps,results] =  perform_rest_regression(data_path,type,AllJID,rest_subjects)
%        - rest_correlation(..., save_results,0)
%        - rest_correlation(..., save_path,'my_dir/correlation')
%
%  Input(s):
%   - data_path = data_path string path to the backfit files (day 1 and day 2)
%   - type = string 'EO' Eyes open or 'EC' Eyes closed
%   - AllJID = struct containing the JIDs
%   - rest_subjects = order of the subjects in ALLJID 
%
%  Optional Input(s):
%     - options.save_results (Default: 1) = logical save output, 0 not to save
%     - options.save_path (Default: current folder) = string or char array path
%
%  Output(s):
%   - results = struct
%       results{state,:} = microstates
%       results{:,fit} = 1 - 'IRLS', 2 - 'OLS'
%   - params_per_pps = cell with microstates parameters for each subject 
%
% Author: R.M.D. Kock

arguments
    data_path;
    type;
    AllJID;
    rest_subjects;
    options.save_results = 1;
    options.save_path = sprintf('%s/regression',data_path);
    options.parameter = 'coverage'
    options.search_condition = 'EEG2'
end
if options.save_results && ~exist(options.save_path, 'dir')
    mkdir(options.save_path); addpath(genpath(save_path))
end
%% prepare regressors
[params_per_pps] = prepare_parameters(data_path,type,options.search_condition);

[x,idx] = sort(rest_subjects(:,1));
AllJID = AllJID(idx,:);
[x2,idx2] = sort(params_per_pps(:,5));
if ~isequal(x,x2)
    error('subjects are not equal')
end
params_per_pps = params_per_pps(idx2,:);
if strcmp(options.parameter,'coverage')
    parameter = cell2mat(params_per_pps(:,3));
elseif strcmp(options.parameter,'duration')
    parameter = cell2mat(params_per_pps(:,1));
elseif strcmp(options.parameter,'occurence')
    parameter = cell2mat(params_per_pps(:,2));
else
    error('select an existing parameter')
end
% run the regression
% regressor = double(table2array(with_jid(:, {'age', 'gender'})));
methods = {'IRLS','OLS'};
n_boot = 1000;
results = cell(10,2);
n_states = size(parameter,2);
% LIMO regression between all pixels of a single JID and microstates regressor
for fit=1:2
    fitMethod = methods{fit};
    for state=1:n_states
        regressor = [log10(parameter(:,state)  + 3.1463e-12 )];
        regressor = parameter(:,state);
        % regressor = [log10(parameter(:,state)  + 3.1463e-12 ),table2array(AllJID(:, 2))];
        [res.val.mask,res.val.p_vals,res.val.mdl,res.val.A, ...
            res.val.B] = singleDayLIMO(regressor, table2array(AllJID(:, 5)), 'FitMethod', fitMethod, 'nBoot', n_boot);
        res.regressor = regressor; 
        res.predictor = table2array(AllJID(:, 5));
        if options.save_results
            save(sprintf('%s/%s_%s_%d_%s.mat',options.save_path,options.parameter,type,state,fitMethod), 'res')
        end
        results{state,fit} = res;
    end
end
if options.save_results
    save(sprintf('%s/params_per_pps_%s.mat',options.save_path,type), 'params_per_pps')
    save(sprintf('%s/%s_results_%s.mat',options.save_path,options.parameter,type), 'results')
end
end