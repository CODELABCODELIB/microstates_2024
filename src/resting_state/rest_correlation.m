function [params_per_pps,Duration,Occurence,Coverage] = rest_correlation(params_per_pps,data_path,type,options)
%% Calculate correlation between resting states of different days 
%
% **Usage:** [params_per_pps,Duration,Occurence,Coverage] = rest_correlation(data_path,type)
%        - rest_correlation(..., save_results,0)
%        - rest_correlation(..., save_path,'my_dir/correlation')
%
%  Input(s):
%   - data_path = data_path string path to the backfit files (day 1 and day 2)
%   - type = string 'EO' Eyes open or 'EC' Eyes closed
%
%  Optional Input(s):
%     - options.save_results (Default: 1) = logical save output, 0 not to save
%     - options.save_path (Default: current folder) = string or char array path
%
%  Output(s):
%   - params_per_pps = cell with microstates parameters for each subject 
%   - Duration = pearson R between file 1 and file 2 for the backfit duration parameter
%   - Occurence = pearson R between file 1 and file 2 for the backfit Occurence parameter
%   - Coverage = pearson R between file 1 and file 2 for the backfit Coverage parameter
%
% Author: R.M.D. Kock

arguments
    params_per_pps;
    data_path;
    type;
    options.save_results logical = 0;
    options.save_path = './';
end
%%
% [params_per_pps] = prepare_parameters(data_path,type,'EEG2');
selected_participants = params_per_pps(:,5);
single_files = unique(selected_participants);
n_subjects = length(single_files);

Duration = nan(n_subjects,1);
Occurence = nan(n_subjects,1);
Coverage = nan(n_subjects,1);
for pp=1:n_subjects
    tmp = find(contains(selected_participants, single_files{pp}));
    if length(tmp) == 2
        P = corrcoef(cell2mat(params_per_pps(tmp(1),1)),cell2mat(params_per_pps(tmp(2),1)));
        %     Duration(pp) = P(2)^2;
        Duration(pp) = P(2);
        P = corrcoef(cell2mat(params_per_pps(tmp(1),2)),cell2mat(params_per_pps(tmp(2),2)));
        %     Occurence(pp) = P(2)^2;
        Occurence(pp) = P(2);
        P = corrcoef(cell2mat(params_per_pps(tmp(1),3)),cell2mat(params_per_pps(tmp(2),3)));
        %     Coverage(pp) = P(2)^2;
        Coverage(pp) = P(2);
    elseif length(tmp) == 3
        P = corrcoef(unique(cell2mat(params_per_pps(tmp,1)), 'rows', 'first'));
        Duration(pp) = P(2);
        P = corrcoef(unique(cell2mat(params_per_pps(tmp,2)), 'rows', 'first'));
        %     Occurence(pp) = P(2)^2;
        Occurence(pp) = P(2);
        P = corrcoef(unique(cell2mat(params_per_pps(tmp,3)), 'rows', 'first'));
        %     Coverage(pp) = P(2)^2;
        Coverage(pp) = P(2);
    end
end
Duration = Duration(~isnan(Duration));
Occurence = Occurence(~isnan(Occurence));
Coverage = Coverage(~isnan(Coverage));
if options.save_results
    save(sprintf('%s/params_per_pps_%s.mat',options.save_path,type), 'params_per_pps')
    save(sprintf('%s/Duration_%s.mat',options.save_path,type), 'Duration')
    save(sprintf('%s/Occurence_%s.mat',options.save_path,type), 'Occurence')
    save(sprintf('%s/Coverage_%s.mat',options.save_path,type), 'Coverage')
end
end