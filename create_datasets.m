function create_datasets( row_id, col_id, site_id, num_frames, if_ERK, filename )

% define parameters
all_CDK2_traces = nan(0, num_frames); all_H2B_traces = all_CDK2_traces; all_other_traces = all_CDK2_traces;

% load traces
for i_row = 1:length(row_id)
    for i_col = 1:length(col_id)
        for i_site = 1:length(site_id)
            % load datasets
            try
                h = load(['../results/', num2str(row_id(i_row)), '_', num2str(col_id(i_col)), '_', num2str(site_id(i_site)), '/signals.mat']);
                all_signals = h.all_signals{row_id(i_row), col_id(i_col), site_id(i_site)};
            catch
                continue;
            end

            % genealogy
            genealogy = nan(length(all_signals), 1);
            for i=1:length(all_signals)
                genealogy(cell2mat(all_signals{i}.daughters)) = i;
            end

            % aggregate traces
            for i=1:length(all_signals)
                if isnan(all_signals{i}.ellipse_id(end))
                    continue;
                end
                curr_H2B_traces = nan(1, num_frames); curr_CDK2_traces = curr_H2B_traces; 
                curr_other_traces = curr_H2B_traces; curr_id = i;
                while 1
                    first_id = find(~isnan(all_signals{curr_id}.ellipse_id), 1, 'first');
                    last_id = find(~isnan(all_signals{curr_id}.ellipse_id), 1, 'last');
                    curr_H2B_traces(first_id:last_id) = all_signals{curr_id}.H2B_nuc_mean(first_id:last_id);
                    curr_CDK2_traces(first_id:last_id) = all_signals{curr_id}.CDK2_cytoring_mean(first_id:last_id) ./ all_signals{curr_id}.CDK2_nuc_mean(first_id:last_id);
                    if (if_ERK)
                        curr_other_traces(first_id:last_id) = all_signals{curr_id}.ERK_cytoring_mean(first_id:last_id) ./ all_signals{curr_id}.ERK_nuc_mean(first_id:last_id);
                    else
                        curr_other_traces(first_id:last_id) = all_signals{curr_id}.p38_cytoring_mean(first_id:last_id) ./ all_signals{curr_id}.p38_nuc_mean(first_id:last_id);
                    end

                    if (first_id == 1)
                        all_H2B_traces = cat(1, all_H2B_traces, smooth(curr_H2B_traces)');
                        all_CDK2_traces = cat(1, all_CDK2_traces, smooth(curr_CDK2_traces)');
                        all_other_traces = cat(1, all_other_traces, smooth(curr_other_traces)');
                        break;
                    end
                    curr_id = genealogy(curr_id);
                    if isnan(curr_id)
                        break;
                    end     
                end
            end
        end
    end
end

% filter traces
low_CDK2_range = [0.2, 0.75]; high_CDK2_range = [1.25, 3]; max_1_frac = 0.4;
valid_cell_id = find(min(all_CDK2_traces, [], 2) >= low_CDK2_range(1) & min(all_CDK2_traces, [], 2) <= low_CDK2_range(2) & ...
    max(all_CDK2_traces, [], 2) >= high_CDK2_range(1) & max(all_CDK2_traces, [], 2) <= high_CDK2_range(2) & ...
    mean(all_CDK2_traces > 0.9 & all_CDK2_traces < 1.1, 2) <= max_1_frac);
all_CDK2_traces = all_CDK2_traces(valid_cell_id, :);
all_H2B_traces = all_H2B_traces(valid_cell_id, :);
all_other_traces = all_other_traces(valid_cell_id, :);

% save datasets
save(['dataset/', filename, '.mat'], 'all_CDK2_traces', 'all_H2B_traces', 'all_other_traces');

end