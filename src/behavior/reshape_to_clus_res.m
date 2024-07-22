function [reshaped_data] = reshape_to_clus_res(cat_params,data)
%% Reshape a dataset to the clustering results shape 
%
% **Usage:**
%   - [all_data] = reshape_to_clus_res(ALLEEG,data)
%
% Input(s):
%    ALLEEG = EEG structs of all participants with microstate field
%    data = data to be reshaped (shape : number of participants x 1)
%
% Output(s):
%   reshaped_data = reshaped input dataset
%
% Author: R.M.D. Kock, Leiden University, 04/12/2023

arguments
    cat_params;
    data cell = {};
end
%%
% [jid_state_params] = create_jid_states(ALLEEG);
reshaped_data = [];
for pp =1:length(cat_params.jid_state_params)
    has_jid = cat_params.jid_state_params{pp,5};
    if ~isempty(data)
        reshaped_data = cat(2,reshaped_data,data{pp}(has_jid));
    end
end
end