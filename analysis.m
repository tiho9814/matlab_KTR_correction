%% preparation
% parameters
all_row_id = {1:8, 1:8, 1:8, 1:8, 1:8, 1:8, 1:8, 1:8, 1:8, 1:8, 1:8, 1:8};
all_col_id = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
all_ref_frames = {64, 64, [64, 96], [64, 96], [64, 96], [64, 96]};
all_filenames = {'Dab', 'Dab_corrected', 'Dab_Tra', 'Dab_Tra_corrected', 'Dab_PF3600', 'Dab_PF3600_corrected'};
%all_ref_frames = {64, 64};
%all_filenames = {'Dab', 'Dab_corrected'};
if_ERK = ones(1, length(all_row_id));
mkdir('dataset'); mkdir('heatmap');

% % create datasets
% for i=1:length(all_row_id)
%     create_datasets( all_row_id{i}, all_col_id{i}, 1:2, 192, if_ERK(i), all_filenames{i} );
% end

% heatmap
for i=1:length(all_row_id)
    plot_heatmap( 4, if_ERK(i), all_ref_frames{i}, all_filenames{i} );
end
