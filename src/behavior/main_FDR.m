function [stats,stats_unadjusted] = main_FDR(data_path, options)
arguments
    data_path;
    options.save_res logical = 1;
    options.save_path = '.';
end
%% load datasets
folder_contents = dir(data_path);
files = {folder_contents.name};
all_files = files(cellfun(@(x) contains(x, 'jid_state_boots'),files));
%% preset variables
stats = cell(length(all_files), 10);
stats_unadjusted = cell(length(all_files), 10);
%% perform FDR
for pp=1:length(all_files)
    tmp_file_name = split(all_files{pp}, 'boots',1);
    % load original jid_state
    load(sprintf('%s/%s',data_path,strcat(tmp_file_name{1},tmp_file_name{2}(2:end))))
    % load jid state boot
    load(sprintf('%s/%s',data_path,all_files{pp}))
    for state=1:length(jid_state_original)
        [stats{pp,state}, stats_unadjusted{pp,state}] = FDR(jid_state_boot, jid_state_original, state);
    end 
end
if options.save_res
    save(sprintf('%s/stats.mat',options.save_path), 'stats')
    save(sprintf('%s/stats_unadjusted.mat',options.save_path), 'stats_unadjusted')
end
end
