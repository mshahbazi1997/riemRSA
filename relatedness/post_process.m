clear; clc; close all
%% === constants ===
voxels = 10:2:50;

sub = 1;
save_addr = 'Result1';
ROIs = [1:10];
n_ROIs = length(ROIs);

%% === load data
meanAcc.RSM.spea = zeros(length(voxels), n_ROIs);
meanAcc.RSM.pear = zeros(length(voxels), n_ROIs);
meanAcc.RSM.kend = zeros(length(voxels), n_ROIs);
meanAcc.RSM.riem = zeros(length(voxels), n_ROIs);
meanAcc.RSM.frob = zeros(length(voxels), n_ROIs);
meanAcc.G.spea = zeros(length(voxels), n_ROIs);
meanAcc.G.pear = zeros(length(voxels), n_ROIs);
meanAcc.G.kend = zeros(length(voxels), n_ROIs);
meanAcc.G.riem = zeros(length(voxels), n_ROIs);
meanAcc.G.frob = zeros(length(voxels), n_ROIs);
meanAcc.G.CKA = zeros(length(voxels), n_ROIs);

count = 0;
for ROI = ROIs
    count = count + 1;
    fileName = sprintf('relatedness_sub0%d_ROI%d.mat', sub, ROI);
    load([save_addr, filesep, fileName], 'acc')
    
    meanAcc.RSM.spea(:, count) = smooth(mean(acc.RSM.spea))';
    meanAcc.RSM.pear(:, count) = smooth(mean(acc.RSM.pear))';
    meanAcc.RSM.kend(:, count) = smooth(mean(acc.RSM.kend))';
    meanAcc.RSM.riem(:, count) = smooth(mean(acc.RSM.riem))';
    meanAcc.RSM.frob(:, count) = smooth(mean(acc.RSM.frob))';
    meanAcc.G.spea(:, count) = smooth(mean(acc.G.spea))';
    meanAcc.G.pear(:, count) = smooth(mean(acc.G.pear))';
    meanAcc.G.kend(:, count) = smooth(mean(acc.G.kend))';
    meanAcc.G.riem(:, count) = smooth(mean(acc.G.riem))';
    meanAcc.G.frob(:, count) = smooth(mean(acc.G.frob))';
    meanAcc.G.CKA(:, count) = smooth(mean(acc.G.CKA))';
    
end

%% === plot ===
% == RSM ==
Fig1 = figure(1);
hold on

plot_woth_conf(voxels, meanAcc.RSM.riem, 'b', n_ROIs)
plot_woth_conf(voxels, meanAcc.RSM.frob, 'k', n_ROIs)
plot_woth_conf(voxels, meanAcc.RSM.spea, 'r', n_ROIs)
plot_woth_conf(voxels, meanAcc.RSM.pear, 'g', n_ROIs)
plot_woth_conf(voxels, meanAcc.RSM.kend, 'm', n_ROIs)


xlabel('# Response channels', 'FontSize', 12); ylabel('Proportion of significant cases', 'FontSize', 12)
legend('Riemannian', 'Frobenius', 'Spearman', 'Pearson', 'Kendall''s \tau', 'FontSize', 12,'Location','southeast')
title('RSM', 'FontSize', 14, 'FontWeight', 'bold')
% set(gcf,'Position',[100 100 300 200])



print(Fig1,[save_addr, filesep, 'RSM_relatedness'],'-dpdf','-r1000')
% saveas(gcf, [save_addr, filesep, 'RSM_relatedness.svg'])
fprintf('Saved in: %s\n', save_addr)

% == G ==

Fig2 = figure(2);
hold on

plot_woth_conf(voxels, meanAcc.G.riem, 'b', n_ROIs)
plot_woth_conf(voxels, meanAcc.G.frob, 'k', n_ROIs)
plot_woth_conf(voxels, meanAcc.G.spea, 'r', n_ROIs)
plot_woth_conf(voxels, meanAcc.G.pear, 'g', n_ROIs)
plot_woth_conf(voxels, meanAcc.G.kend, 'm', n_ROIs)
plot_woth_conf(voxels, meanAcc.G.CKA, 'c', n_ROIs)

xlabel('# Response channels', 'FontSize', 12); ylabel('Proportion of significant cases', 'FontSize', 12)
legend('Riemannian', 'Frobenius', 'Spearman', 'Pearson', 'Kendall''s \tau', 'CKA', 'FontSize', 12, 'Location','southeast')
title('2^{nd}-moment', 'FontSize', 14, 'FontWeight', 'bold')
%set(gcf,'Position',[100 100 300 200])


print(Fig2,[save_addr, filesep, 'G_relatedness'],'-dpdf','-r1000')
% saveas(gcf, [save_addr, filesep, 'G_relatedness.svg'])
fprintf('Saved in: %s\n', save_addr)

% %% === newer plot ===
% 
% 
% imagesc(flip(meanAcc.RSM.frob',2))
% axis square;
% colormap default
% caxis([0 1])
% xticks([1, 21]);
% xticklabels([voxels(end), voxels(1)])
% yticks([1, 10])
% yticklabels([1, 10])
% 
% 
% xlabel('voxel num.', 'FontSize', 14)
% ylabel('region num.', 'FontSize', 14)

%% === plot function ===

function plot_woth_conf(x_axis, meter, color, n_sample)
    meter_mean = mean(meter, 2)';
    not_nan = ~isnan(meter_mean);

    plot(x_axis(not_nan), meter_mean(not_nan), color, 'LineWidth', 2)
    
    meter_conf = std(meter, [], 2)'/sqrt(n_sample);
        
    plt1 = fill([x_axis(not_nan), fliplr(x_axis(not_nan))], ...
        [meter_mean(not_nan) + meter_conf(not_nan), fliplr(meter_mean(not_nan) - meter_conf(not_nan))], ...
        color, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    alpha(plt1, 0.2)
end