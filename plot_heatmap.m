function plot_heatmap(frame_rate, if_ERK, ref_frame, filename)

% load data
h = load(['dataset/', filename, '.mat']);
all_CDK2_traces = h.all_CDK2_traces; all_other_traces = h.all_other_traces;

% sort by similarity of CDK2 traces
D = pdist(all_CDK2_traces, 'spearman');
leafOrder = optimalleaforder(linkage(D, 'average'), D);
all_CDK2_traces = all_CDK2_traces(leafOrder, :);
all_other_traces = all_other_traces(leafOrder, :);

% colormap
temp = redbluecmap;
newCmap = imresize(temp, [64, 3]); newCmap = min(max(newCmap, 0), 1);
newCmap_cropped = imresize(temp(3:end, :), [64, 3]); newCmap_cropped = min(max(newCmap_cropped, 0), 1);

% plot data, CDK2
h = figure(1); T = (1:size(all_CDK2_traces, 2))/4;
imagesc(T, 1:size(all_CDK2_traces, 1), all_CDK2_traces); hold on; 
xlim([0, max(T)]); ylim([1, size(all_CDK2_traces, 1)]);
xlabel('Time (Hours)'); xticks(0:12:48); set(gca, 'ydir', 'reverse');
colormap(newCmap_cropped); caxis([0.25, 1.75]); 
cb = colorbar; set(cb, 'Ticks', 0.25:0.25:1.75);
ylabel('CDK2 Activity');
for i=1:length(ref_frame)
    plot([1, 1]*ref_frame(i)/frame_rate, [1, size(all_CDK2_traces, 1)], 'k', 'linewidth', 2); 
end
h.Renderer = 'Painters'; h.PaperUnits = 'inches';
h.PaperPosition = [0, 0, 4, 4]; h.PaperSize = [4, 4];
print(h, '-dpng', '-r600', ['heatmap/', filename, '_CDK2.png']);
close(h);

% plot data, sensor 2
h = figure(1); 
imagesc(T, 1:size(all_other_traces, 1), all_other_traces); hold on; 
xlim([0, max(T)]); ylim([1, size(all_other_traces, 1)]);
xlabel('Time (Hours)'); xticks(0:12:48); set(gca, 'ydir', 'reverse');
%colormap(newCmap); cb = colorbar; 
colormap(jet); cb = colorbar; 
if (if_ERK)
    %caxis([0.5, 1.5]); set(cb, 'Ticks', 0.5:0.25:1.5);
    caxis([0.2, 2]); set(cb, 'Ticks', 0.2:0.2:2);
    ylabel('ERK Activity');
else
    caxis([0.25, 0.75]); set(cb, 'Ticks', 0.25:0.1:0.75);
    ylabel('p38 Activity');
end

for i=1:length(ref_frame)
    plot([1, 1]*ref_frame(i)/frame_rate, [1, size(all_other_traces, 1)], 'k', 'linewidth', 2); 
end
h.Renderer = 'Painters'; h.PaperUnits = 'inches';
h.PaperPosition = [0, 0, 4, 4]; h.PaperSize = [4, 4];
print(h, '-dpng', '-r600', ['heatmap/', filename, '_other.png']);
close(h);

end