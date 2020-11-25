clear
clc
close all

%%
tit = 'F1-Score V.S. Metrics';

load_addr = 'post_results';

n_sub = 20;

measures = {'Riem', 'Frob', 'Spea', 'Pear', 'Kend'};

measures_complete = { 'Riemannian', 'Frobenius', 'Spearman', 'Pearson', 'Kendall''s \tau'};

measure_colors = {'b', 'k', 'r', 'g', 'm'};

n_seed = 8;

%%
n_meas = size(measures, 2);

% init.
F1_seed_meas = nan(n_seed, n_meas);
p = nan(n_meas, n_meas);

F1_cell = load([load_addr, filesep, sprintf('F1_%d.mat', n_sub)]).F1_all;

%%
for i_meas = 1 : n_meas
    measure_i = measures{i_meas};
    
    F1_seed_meas(:, i_meas) = F1_cell.(measure_i);

    for j_meas = 1 : n_meas
        measure_j = measures{j_meas};

        p(i_meas, j_meas) = signrank(F1_cell.(measure_i), F1_cell.(measure_j), 'tail', 'right')/2;
    end
end
    

%%
fig = figure
hold on

h_box = boxplot(F1_seed_meas, 'symbol', '');
set(h_box, {'linew'}, {2})

%%
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');

ccc = 0.02;
riem_idx = 1;

max_yt = max(F1_seed_meas(:));
min_yt = min(F1_seed_meas(:));

index = 1;


for j_measure = 1 : n_meas
    if j_measure == riem_idx
        continue
    end
    
    if p(riem_idx, j_measure) > 0.05
        continue
    end
    
    
    line_y = [1, 1]*max_yt*(1 + (index*ccc)/max_yt);

    line_x = xt([riem_idx, j_measure]); 
    
    plot(line_x, line_y, '-k', 'LineWidth', 1.5)
    

    index = index + 1;
    
end

ylim([min_yt*(1 - ccc/min_yt), max_yt*(1 + (index*ccc)/max_yt)])

ylabel('F1-score')

set(gca, 'xtick', 1 : n_meas, 'xticklabel', measures_complete, ...
    'XTickLabelRotation', -45, ...
    'FontSize', 12, ...
    'TickLabelInterpreter', 'tex')

%%
print(fig, fullfile([load_addr, filesep,' F1.pdf']), '-dpdf', '-r1000')