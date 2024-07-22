function [labels_renamed] = rename_microstates_labels(labels,options)
arguments
    labels
    options.dictionary_labels = [6,5,1,8,10,12,11,2,3,9,4,7]; 
    options.microstates_labels={'A','B','C','D','E', 'F', 'G', 'H', 'I', 'J', 'K', 'L'};
end
labels_renamed = cell(length(labels),1);
for i=1:length(labels)
    current_label = labels(i);    
    labels_renamed{i} = options.microstates_labels{options.dictionary_labels == current_label};
end

for i=1:length(labels)
    current_label = labels(i); 
    duplicate_indices = find(labels == current_label);
    if length(duplicate_indices) > 1 && numel(char(labels_renamed{i})) <= 1
        for d=1:length(duplicate_indices)
            labels_renamed{duplicate_indices(d)} = strcat(labels_renamed{duplicate_indices(d)}, int2str(d));
        end
    end
end
end