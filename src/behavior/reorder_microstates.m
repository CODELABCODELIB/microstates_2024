function [indices,reordered_microstates_labels] = reorder_microstates(individual_microstates_labels,options)
%% Reorder the microstates at individual level to population level microstates
%
% **Usage:**
%   [indices,reordered_microstates_labels] = reorder_microstates(individual_microstates_labels)
%                   reorder_microstates(...,[1 2 3 4 5])
%
% Input(s):
%    - data_array = current microstates labels
%
% Optional input parameter(s):
%   - dictionary_labels = order of microstates as in paper
%
% % Output(s):
%   - indices = the indices of the reordered labels
%   - reordered_microstates_labels = reordered labels
%
% Author(s)
% Ruchella Kock, Leiden University, 15/11/2023

arguments
    individual_microstates_labels
    options.dictionary_labels = [6,5,1,8,10,12,11,2,3,9,4,7]; 
end
% copying the data array so that we have our original with us  
data_array_copy=individual_microstates_labels; 
n=numel(options.dictionary_labels); 
% array to store the indices of the nearest values 
indices = nan(1,length(individual_microstates_labels));
count=1;
for k = 1:n 
    % abs(data_array_copy - test_array(k)) will give the difference array 
    % and in that difference array we find the minimum index 
    % Note we need the index because we also have to mark the element in 
    % data_array as NaN 
    differences = abs(data_array_copy-options.dictionary_labels(k));
%     [smallest_diff]=min(differences); 
    tmp_ind = find(differences==0);
    if isempty(tmp_ind)
        continue
    end
    indices(count:count+length(tmp_ind)-1) = tmp_ind;
    data_array_copy(tmp_ind)=NaN; 
    count = count+length(tmp_ind);
end 
reordered_microstates_labels = individual_microstates_labels(indices);
end