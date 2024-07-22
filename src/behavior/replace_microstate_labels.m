function [microstates] = replace_microstate_labels(microstates,reordered_microstates_labels)
unique_states = unique(microstates);
logical_indices = zeros(length(unique_states),length(microstates));
% first get the indices do not replace yet because this may lead to 
% replacing of replaced values
for i=1:length(unique_states)
    logical_indices(i,:) = microstates == unique_states(i);
end
logical_indices = logical(logical_indices);
for i=1:length(unique_states)
    microstates(logical_indices(i,:)) = reordered_microstates_labels(i);
end
