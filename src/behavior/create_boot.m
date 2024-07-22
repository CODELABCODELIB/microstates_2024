function [boot_data] = create_boot(labels,state)
shuffled_labels = labels(randperm(length(labels)));
shuffled_labels_indx = find(shuffled_labels == state);
boot_data = datasample(shuffled_labels_indx,length(shuffled_labels_indx));
end