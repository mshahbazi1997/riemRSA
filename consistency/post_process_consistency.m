clear; clc; close all;

%% === Settings ===
use_rsm = false;

load_addr = 'results/Haxby';
save_addr = [load_addr, filesep, 'post_process'];
mkdir(save_addr)
file_name = ['bias_corrected_', int2str(use_rsm)];

roi_labels = {'V1(L)', 'V1(R)', 'V2(L)', 'V2(R)', 'V3(L)', 'V3(R)', 'V4(L)', 'V4(R)', 'V5(L)', 'V5(R)'};
measures = {'Spearman', 'Pearson', 'Riemannian', 'Frobenius', 'Kendall''s \tau', 'CKA'};
measures_order_all = [3 1 2 5 4 6];

%% === Load results ===
load([load_addr, filesep, file_name])
% bias_corrected: [roi x comb_run x sub x measure]

%% === Process settings ===
[n_roi, ~, n_sub, n_measure] = size(bias_corrected);

measures_order = measures_order_all(1:n_measure);

bias_corrected = bias_corrected(:, :, :, measures_order);

%% === Consistency score ===
consistency_roi_measure = nan(n_roi, n_measure);

for i_measure = 1 : n_measure
    bias_corrected_roi_sub = nan(n_roi, n_sub);    
    for i_sub = 1 : n_sub
        bias_corrected_roi_sub(:, i_sub) = mean(squeeze(bias_corrected(:, :, i_sub, i_measure)), 2);
    end
    bias_corrected_roi = mean(bias_corrected_roi_sub, 2)./std(bias_corrected_roi_sub, [], 2)/sqrt(n_sub);
    
    consistency_roi_measure(:, i_measure) = bias_corrected_roi;
end

%%
% figure
% hold on
% 
% roi_dc = 0.2;
% label_ys = [];
% for roi = 1 : n_roi
%     if roi == 1
%         dc_i = 0;
%     else
%         dc_i = dc_i + max(consistency_roi_measure(roi - 1, :)) + roi_dc;
%     end
%     
%     plot(1 : n_measure, dc_i + consistency_roi_measure(roi, :), 'LineWidth', 3, 'Color', 'k')
%     
%     plot(1 : n_measure, dc_i*ones(1, n_measure), 'r--')
%     
%     label_ys = [label_ys, dc_i + consistency_roi_measure(roi, 1)];
% end
% set(gcf, 'position',[0, 0, (n_measure - 1)*68, 450])
% set(gca, 'ytick', label_ys, 'yticklabel', roi_labels, 'YTickLabelRotation', 0, 'FontSize', 14)
% set(gca, 'xtick', 1 : n_measure, 'xticklabel', measures(measures_order), 'xTickLabelRotation', -45, 'FontSize', 14)
% 
% if use_rsm
%     title('RSM', 'FontSize', 16)
% else
%     title('G', 'FontSize', 16)
% end

% saveas(gcf, [save_addr, filesep, 'figure_roi_measure_', int2str(use_rsm)], 'svg')

%%
pval_matrix = nan(n_measure, n_measure);
for i_measure = 1 : n_measure
    for j_measure = 1 : n_measure
        [pval, ~, ~] = signrank(consistency_roi_measure(:,i_measure), consistency_roi_measure(:, j_measure), 'tail', 'right');
        pval_matrix(i_measure, j_measure) = pval;
    end
end

%%
d = consistency_roi_measure; %randi(10, 20, 3);


pThreshold = FDRthreshold(pval_matrix(1,2:end),0.05,1);

Fig = figure;

x = repmat((1:n_measure)',1, n_roi); 

plot(x,d','-o', 'color',[0.9 0.9 0.9],'MarkerEdgeColor',[0.6 0.6 0.6], 'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',2)

hold on


h_box = boxplot(d, 'symbol', '');
set(h_box, {'linew'}, {2})

%
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');

hold on
ccc = 0.30;
riem_idx = 1;
% astrisk_dc = 0.35;
astrisk_txt = 0.04;
% astrisk_step = 0.1;

max_yt = prctile(d(:,1),75);%max(yt);

index = 0;


for i_measure = 1 : n_measure
    if i_measure == riem_idx
        continue
    end
    
    if pval_matrix(riem_idx, i_measure)>pThreshold
        continue
    end
    
    
    line_y = [1, 1]*max_yt*(1 + (index*ccc)/max_yt);
%     text_y = max_yt*(1 + (astrisk_txt + index*ccc)/max_yt);
    line_x = xt([riem_idx, i_measure]); 
%     text_x = mean(xt([riem_idx, i_measure]));
    
%     if i_measure > riem_idx
%         line_y = line_y - 0.1*astrisk_step;%max_yt*astrisk_step;
%         text_y = text_y - 0.1*astrisk_step;%max_yt*astrisk_step;
%     end
    
    
    plot(line_x, line_y, '-k', 'LineWidth', 1.5)
    
%     text(text_x, text_y, sprintf('*%.3f', pval_matrix(riem_idx, i_measure)), 'FontSize', 13);
%     text(text_x, text_y, '*', 'FontSize', 13);
    
    index = index + 1;
end

hold off

if use_rsm
%     title('RSM', 'FontSize', 14)
    
    ylim([0, 2.6])
else
%     title('2^{nd}-moment', 'FontSize', 14)
    ylim([0, 2.6])
end


set(gca, 'xtick', 1 : n_measure, 'xticklabel', measures(measures_order), 'XTickLabelRotation', -45, 'FontSize', 12 ...
    , 'TickLabelInterpreter', 'tex')

% set(gcf, 'position',[0, 0, (n_measure - 1)*68, 320])

% ylim([0, 1.9])

% saveas(gcf, [save_addr, filesep, 'figure_boxplot_', int2str(use_rsm)], 'svg')%'epsc'
print(Fig,[save_addr, filesep, 'figure_boxplot_', int2str(use_rsm)],'-dpdf','-r1000')

% print(gcf, [save_addr, filesep, 'figure_boxplot_', int2str(use_rsm), '.png'],'-dpng','-r2000')

% %%
% figure 
% 
% x = kron((1:n_measure)',ones(n_roi, 1));
% 
% % violinplot(x, consistency_roi_measure, 'edgecolor', 'none', 'facecolor', [0.5, 0.5, 0.5], 'meancolor', 'r');
% myboxplot(x, consistency_roi_measure);
% 
% %%
% ee = measures(measures_order);
% 
% ax = gca;
% ax.FontSize = 14; ax.LabelFontSizeMultiplier = 1; ax.TitleFontSizeMultiplier = 12;
% xtickangle(ax,45)
% for l = 1:numel(ax.XTickLabel)
%     ax.XTickLabel{l} = ee{l};
% end