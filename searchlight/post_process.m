clear
clc
close all

%% ===== Constants =====
n_sub_in_split = 20;
n_sub_tot = 160;

%% ===== Addresses =====
curr_addr = pwd; 

depend_path = 'dependencies';
         
load_addr = 'results';
save_addr = ['post_', load_addr];
mkdir(save_addr)

toolbox_path = 'rsatoolbox-selected';
addpath(genpath(toolbox_path));

low_level_funcs_path = '../low_level_functions';
addpath(low_level_funcs_path)

%% ===== Load anatomical volume and its mask =====
load([load_addr, filesep, 'Mask.mat'])
load([depend_path, filesep, 'sampleMask_org.mat'])
load([depend_path, filesep, 'anatomy.mat'])

%% ===== Load maps =====
for measure_cell = {'Riem', 'Frob'}
    measure = measure_cell{1};
    [p_maps.(measure), bc_maps.(measure)] = get_p_map_ttest( ...
        [load_addr, filesep, 'maps'], 1:n_sub_tot, 'dist', measure);
end

for measure_cell = {'Spea', 'Pear', 'Kend'}
    measure = measure_cell{1};
    [p_maps.(measure), bc_maps.(measure)] = get_p_map_ttest( ...
        [load_addr, filesep, 'maps'], 1:n_sub_tot, 'corr', measure);
end

%% ===== Brain demo =====
fig = figure;

th_new = BonferroniThreshold(p_maps.Riem, 0.05);

y_th = zeros(size(p_maps.Riem));
y_th(p_maps.Riem <= th_new) = 1;

brainVol = addRoiToVol(map2vol(anatVol), mask2roi(m), [0.5 0.5 1], 2);
brainVol_y = addBinaryMapToVol(brainVol, y_th.*m, [1 1 0.3]);
brainVol_sr = addRoiToVol(brainVol_y, mask2roi(Mask.*m), [1 0.35 0.35], 2);
showVol(brainVol_sr, '', 1)

%%
print(fig, fullfile([save_addr, filesep,' brain_demo.pdf']), '-dpdf', '-r1000')

%% ===== Test on splits ====
n_split = floor(size(bc_maps.Riem, 1)/n_sub_in_split);
th = 0.05;

% init
for measure_cell = {'Riem', 'Frob', 'Spea', 'Pear', 'Kend'}
    measure = measure_cell{1};
    F1_all.(measure) = nan(1, n_split);
    percision_all.(measure) = nan(1, n_split);
    recall_all.(measure) = nan(1, n_split);
end

%%
all_p_maps = {};

for i_split = 1 : n_split
    disp(i_split)
    
    %% ===== T-Test =====    
    subs = 1 + (i_split - 1)*n_sub_in_split : i_split*n_sub_in_split;
    
    for measure_cell = {'Riem', 'Frob', 'Spea', 'Pear', 'Kend'}
        tic 
        
        measure = measure_cell{1};
        
        [~, p_map_ttest] = ttest(bc_maps.(measure)(subs, :, :, :), 0, 'Tail', 'right');
        all_p_maps.(measure) = squeeze(p_map_ttest);

        toc
    end

    %% ===== Thresholding =====
    for measure_cell = {'Spea', 'Pear', 'Kend', 'Riem', 'Frob'}
        measure = measure_cell{1};

        th_new = FDRthreshold(all_p_maps.(measure), th);

        y_pr = zeros(size(all_p_maps.(measure)));
        y_pr(all_p_maps.(measure) <= th_new) = 1;

        TP.(measure) = sum(Mask(:).*y_pr(:));
        TN.(measure) = sum(~Mask(:).*~y_pr(:));
        FP.(measure) = sum(~Mask(:).*y_pr(:));
        FN.(measure) = sum(Mask(:).*~y_pr(:));

        accuracy.(measure) = (TN.(measure) + TP.(measure))/(TN.(measure) + TP.(measure) + FP.(measure) + FN.(measure));
        percision.(measure) = TP.(measure)/(TP.(measure) + FP.(measure));
        recall.(measure) = TP.(measure)/(TP.(measure) + FN.(measure));
        F1.(measure) = 2*(recall.(measure).*percision.(measure))/(recall.(measure) + percision.(measure));

        F1_all.(measure)(i_split) = F1.(measure);
        percision_all.(measure)(i_split) = percision.(measure);
        recall_all.(measure)(i_split) = recall.(measure);
    end

end

%%
figure
hold on

plot_with_error_bar(1, F1_all.Riem, 'b')
plot_with_error_bar(2, F1_all.Frob, 'k')
plot_with_error_bar(3, F1_all.Spea, 'r')
plot_with_error_bar(4, F1_all.Pear, 'g')
plot_with_error_bar(5, F1_all.Kend, 'm')


xlim([0.5, 7])

set(gca, 'XTick', 1 : 5, 'XTickLabel', {'Riem', 'Frob', 'Spea', 'Pear', 'Kend'}, 'xTickLabelRotation', -45, 'FontSize', 18)

xlabel('p-value threshold', 'FontSize', 16)

title('F1', 'FontSize', 18)

legend('Riemannian', 'Frobenius', 'Spearman', 'Pearson', 'Kendall''s tau', 'FontSize', 12)

%%
save([save_addr, filesep, sprintf('F1_%d.mat', n_sub_in_split)], 'F1_all')

%%
function plot_with_error_bar(x_axis, meter, color)
    meter_mean = mean(meter(~isnan(meter)));
    
    meter_conf = std(meter(~isnan(meter)));
    
    errorbar(x_axis, meter_mean, meter_conf, 's', ...
        'Color', color, ...
        'LineWidth', 2.5, ...
        'CapSize', 24, ...
        'MarkerSize', 14, ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color)
end