clear; clc; close all;

%% === Settings ===
use_rsm = true;

load_addr = 'results/Haxby';
save_addr = [load_addr, filesep, 'post_process'];
mkdir(save_addr)
file_name = ['bias_corrected_', int2str(use_rsm)];

roi_labels = {'', '', '', '', '', '', '', '', '', ''};
measures = {'Spearman', 'Pearson', 'Riemannian', 'Frobenius', 'Kendall''s \tau', 'CKA'};
measures_order_all = [3 1 2 5 4 6];

%% === Load results ===
load([load_addr, filesep, file_name])
% bias_corrected: [roi x comb_run x sub x measure]

%% === Process settings ===
[n_roi, ~, n_comb_of_runs, n_sub, n_measure] = size(bias_corrected);

measures_order = measures_order_all(1:n_measure);

bias_corrected = bias_corrected(:, :, :, :, measures_order);

%% === Discriminability score ===

bc_roi_roi_sub_measure = squeeze(mean(bias_corrected, 3));

disc_roi_roi_measure = squeeze(mean(bc_roi_roi_sub_measure, 3)./std(bc_roi_roi_sub_measure, [], 3)/sqrt(n_sub));

%%
mask = ones(n_roi);
mask = triu(mask, 1) == 1;

pval_matrix = nan(n_measure);
disc_combroi_measure = nan(sum(mask(:)), n_measure);
for i_measure = 1 : n_measure
    disc_roi_roi_i = disc_roi_roi_measure(:, :, i_measure);
    disc_combroi_measure(:, i_measure) = disc_roi_roi_i(mask);
    for j_measure = 1 : n_measure
        disc_roi_roi_j = disc_roi_roi_measure(:, :, j_measure);
        [pval, ~, ~] = signrank(disc_roi_roi_i(mask), disc_roi_roi_j(mask), 'tail', 'right');
        pval_matrix(i_measure, j_measure) = pval;
    end
end

%%
d = disc_combroi_measure; %randi(10, 20, 3);

Fig = figure;

x = repmat((1:n_measure)',1,size(d,1)); 

plot(x,d','-o', 'color',[0.9 0.9 0.9],'MarkerEdgeColor',[0.6 0.6 0.6], 'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',2)

hold on



h_box = boxplot(d, 'symbol', '');
set(h_box, {'linew'}, {2})



% myboxplot(x, d)
% hold on


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
    
    if pval_matrix(riem_idx, i_measure)>0.05
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
    ylim([0, 3.5])
else
%     title('2^{nd}-moment', 'FontSize', 14)
    ylim([0, 3.5])
end

set(gca, 'xtick', 1 : n_measure, 'xticklabel', measures(measures_order), 'XTickLabelRotation', -45, 'FontSize', 12 ...
    , 'TickLabelInterpreter', 'tex')

% set(gcf, 'position',[0, 0, (n_measure - 1)*68, 320])

% saveas(gcf, [save_addr, filesep, 'figure_boxplot_', int2str(use_rsm)], 'svg')
print(Fig,[save_addr, filesep, 'figure_boxplot_', int2str(use_rsm)],'-dpdf','-r1000')
