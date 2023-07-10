clear vars;
clear all;

%% load Dab data

load('/Users/tim/Desktop/manuscript_ERKKTR_CellSystems/correction_analysis_TH/dataset/Dab.mat');

%% Plot some single-cell traces

% figure(1);
% 
% for j = 1:50
%     subplot(5,10,j);
%       mdl = fitlm(all_CDK2_traces(j,70:192),all_other_traces(j,70:192));
%       plot(mdl, 'Marker', '.');
%       title(num2str(j));
%       xlabel('CDK2 activity');
%       ylabel('ERK activity');
%       %xlim([0.25 1.25]);
%       %ylim([0.25 0.75]);
%       legend off;
% end

figure(2);

for j = 1:10
    subplot(1,10,j);
      plot(all_CDK2_traces(j,:)); hold on; plot(all_other_traces(j,:));
      title(num2str(j));
     % xlabel('Frame');
      ylabel('Repoter activity');
      %xlim([0.25 1.25]);
      ylim([0.25 2]);
      legend off;
end

time = [0.25:0.25:48];

figure(3)

    subplot(1,3,1);
      plot(time, all_CDK2_traces(44,:)); hold on;
      ylim([0.25 2.25]);
      ylabel('CDK2 activity');
      xline(16, '-k')
      xlim([0 48])
      xlabel('Time (h)');
      xticks([0 12 24 36 48])
    yyaxis right;
      plot(time, all_other_traces(44,:));
      ylabel('ERK activity');
      ylim([0.35 1.5])
      
    subplot(1,3,2);
      plot(time, all_CDK2_traces(64,:)); hold on;
      ylim([0.25 2.25]);
      ylabel('CDK2 activity');
      xline(16, '-k')
      xlim([0 48])
      xlabel('Time (h)');
      xticks([0 12 24 36 48])
    yyaxis right;
      plot(time, all_other_traces(64,:));
      ylabel('ERK activity');
      ylim([0.35 1.5])
      
    subplot(1,3,3);
      plot(time, all_CDK2_traces(82,:)); hold on;
      ylim([0.25 2.25]);
      ylabel('CDK2 activity');
      xline(16, '-k')
      xlim([0 48])
      xlabel('Time (h)');
      xticks([0 12 24 36 48])
    yyaxis right;
      plot(time, all_other_traces(82,:));
      ylabel('ERK activity');
      ylim([0.35 1.5])

%% Calculate R squared for a single cell (CDK2 tracedata vs. ERK tracedata) for last 123 frames

CDK2 = all_CDK2_traces(1,70:192);
ERK = all_other_traces(1,70:192);
[b,bint,r,rint,stats] = regress((CDK2)', (ERK)');
R2=stats(1,1);


%% Calculate optimal coefficient for CDK2 vs. (ERK - CDK2*k) for all traces

figure(4);

R2_store = [];

for k = 0:0.001:0.5
    
    for i = 1:6444
       CDK2_k = (all_CDK2_traces(i,70:192)*k);
       ERK = all_other_traces(i,70:192);
       ERK_subtract = (ERK-CDK2_k);
       [b,bint,r,rint,stats] = regress((all_CDK2_traces(i,70:192))',(ERK_subtract)');
       R2_corrected=stats(1,1);
       R2_store(1,i) = R2_corrected;
    end
    
 R2_corrected_mean = mean(R2_store);
    
plot(k,R2_corrected_mean,'.','MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on

yline(0, '-r');
ylim([-2 1]);
xlim([0 0.5]);

end


%% Save corrected dataset file

all_CDK2_traces_k = (all_CDK2_traces * 0.217);

all_other_traces = (all_other_traces - all_CDK2_traces_k);

save(['dataset/Dab_corrected.mat'], 'all_CDK2_traces', 'all_H2B_traces', 'all_other_traces');
