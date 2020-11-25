clear
clc

tic
%% ===== User Inputs =====
low_level_func_path = '../utils/Low-Level Functions';

n_rept = 30;
n_voxels = floor(linspace(10, 100, 50));

do_perm = false;

method = '2nd-moment';

rng(1)

%% ===== Add Path/Load =====
load('dependencies/all_RSMs.mat')

addpath(low_level_func_path)

%% ===== Inline Helper Functions =====
similarity_func = ...
    @(R_a, R_b, pval_perm_func, varargin) out2(pval_perm_func, R_a, R_b, false, varargin);

score_func = @(R_1, R_1_n, R_2, R_2_n, pval_perm_func, varargin) ...
    ((similarity_func(R_1, R_1_n, pval_perm_func, varargin) >= similarity_func(R_1, R_2_n, pval_perm_func, varargin)) + ...
    (similarity_func(R_2, R_2_n, pval_perm_func, varargin) >= similarity_func(R_2, R_1_n, pval_perm_func, varargin)) + ...    
    (similarity_func(R_1, R_1_n, pval_perm_func, varargin) >= similarity_func(R_2, R_1_n, pval_perm_func, varargin)) + ...
    (similarity_func(R_2, R_2_n, pval_perm_func, varargin) >= similarity_func(R_1, R_2_n, pval_perm_func, varargin)))/4;


%% ===== Big Loop =====
h_wb = waitbar(0, 'Processing, Please wait...');

[n_run, n_ROI, n_cond, ~] = size(all_RSMs);

all_scores.Pear = nan(n_run, n_ROI, n_ROI, length(n_voxels));
all_scores.Spea = nan(n_run, n_ROI, n_ROI, length(n_voxels));
all_scores.Kend = nan(n_run, n_ROI, n_ROI, length(n_voxels));
all_scores.Riem = nan(n_run, n_ROI, n_ROI, length(n_voxels));
all_scores.Frob = nan(n_run, n_ROI, n_ROI, length(n_voxels));
all_scores.CKA = nan(n_run, n_ROI, n_ROI, length(n_voxels));

for i_run = 1 : n_run
    waitbar(i_run/n_run, h_wb)
    
    for ROI_1 = 1 : n_ROI
        for ROI_2 = ROI_1 + 1 : n_ROI
            R_1_ideal = squeeze(all_RSMs(i_run, ROI_1, :, :));
            R_2_ideal = squeeze(all_RSMs(i_run, ROI_2, :, :));
            
            for i_voxel = 1 : length(n_voxels)
                n_voxel = n_voxels(i_voxel);
                
                this_run_roi_nvoxel_scores.Pear = nan(1, n_rept);
                this_run_roi_nvoxel_scores.Spea = nan(1, n_rept);
                this_run_roi_nvoxel_scores.Kend = nan(1, n_rept);
                this_run_roi_nvoxel_scores.Riem = nan(1, n_rept);
                this_run_roi_nvoxel_scores.Frob = nan(1, n_rept);
                this_run_roi_nvoxel_scores.CKA = nan(1, n_rept);
                for i_rept = 1 : n_rept
                    B_R_1_1st = mvnrnd(zeros(n_cond, 1), R_1_ideal, n_voxel)';
                    B_R_1_2nd = mvnrnd(zeros(n_cond, 1), R_1_ideal, n_voxel)';
                    
                    B_R_2_1st = mvnrnd(zeros(n_cond, 1), R_2_ideal, n_voxel)';
                    B_R_2_2nd = mvnrnd(zeros(n_cond, 1), R_2_ideal, n_voxel)';
                    
                    % form RSMs
                    R_1_1st = get_representation(B_R_1_1st, method);
                    R_1_2nd = get_representation(B_R_1_2nd, method);
                    R_2_1st = get_representation(B_R_2_1st, method);
                    R_2_2nd = get_representation(B_R_2_2nd, method);
                    
                    % score calculation
                    this_run_roi_nvoxel_scores.Pear(i_rept) = score_func(R_1_1st, R_1_2nd, R_2_1st, R_2_2nd, @pval_perm_Pearson, [], method=='RSM');
                    this_run_roi_nvoxel_scores.Spea(i_rept) = score_func(R_1_1st, R_1_2nd, R_2_1st, R_2_2nd, @pval_perm_Spearman, [], method=='RSM');
                    this_run_roi_nvoxel_scores.Kend(i_rept) = score_func(R_1_1st, R_1_2nd, R_2_1st, R_2_2nd, @pval_perm_Kendall, [], method=='RSM');
                    this_run_roi_nvoxel_scores.Riem(i_rept) = 1 - score_func(R_1_1st, R_1_2nd, R_2_1st, R_2_2nd, @pval_perm_Riemann, [], method=='RSM');
                    this_run_roi_nvoxel_scores.Frob(i_rept) = 1 - score_func(R_1_1st, R_1_2nd, R_2_1st, R_2_2nd, @pval_perm_Frobenius, [], method=='RSM');
                    % cka score calc.
                    Gc_1_1st = get_representation(B_R_1_1st, 'G-centered');
                    Gc_1_2nd = get_representation(B_R_1_2nd, 'G-centered');
                    Gc_2_1st = get_representation(B_R_2_1st, 'G-centered');
                    Gc_2_2nd = get_representation(B_R_2_2nd, 'G-centered');
                    this_run_roi_nvoxel_scores.CKA(i_rept) = score_func(Gc_1_1st, Gc_1_2nd, Gc_2_1st, Gc_2_2nd, @pval_perm_CKA);
                end
                
                all_scores.Pear(i_run, ROI_1, ROI_2, i_voxel) = mean(this_run_roi_nvoxel_scores.Pear);
                all_scores.Spea(i_run, ROI_1, ROI_2, i_voxel) = mean(this_run_roi_nvoxel_scores.Spea);
                all_scores.Riem(i_run, ROI_1, ROI_2, i_voxel) = mean(this_run_roi_nvoxel_scores.Riem);
                all_scores.Frob(i_run, ROI_1, ROI_2, i_voxel) = mean(this_run_roi_nvoxel_scores.Frob);
                all_scores.Kend(i_run, ROI_1, ROI_2, i_voxel) = mean(this_run_roi_nvoxel_scores.Kend);
                all_scores.CKA(i_run, ROI_1, ROI_2, i_voxel) = mean(this_run_roi_nvoxel_scores.CKA);
            end
        end
    end
end

close(h_wb)

toc

%% ===== Overall Score =====
score.Pear = nan(1, length(n_voxels));
score.Spea = nan(1, length(n_voxels));
score.Kend = nan(1, length(n_voxels));
score.Riem = nan(1, length(n_voxels));
score.Frob = nan(1, length(n_voxels));
score.CKA = nan(1, length(n_voxels));

score_std.Pear = nan(1, length(n_voxels));
score_std.Spea = nan(1, length(n_voxels));
score_std.Kend = nan(1, length(n_voxels));
score_std.Riem = nan(1, length(n_voxels));
score_std.Frob = nan(1, length(n_voxels));
score_std.CKA = nan(1, length(n_voxels));

for i_voxel = 1 : length(n_voxels)
    all_scores_nvoxel = all_scores.Pear(:, :, :, i_voxel);
    score.Pear(i_voxel) = mean(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    score_std.Pear(i_voxel) = std(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    
    all_scores_nvoxel = all_scores.Spea(:, :, :, i_voxel);
    score.Spea(i_voxel) = mean(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    score_std.Spea(i_voxel) = std(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    
    all_scores_nvoxel = all_scores.Kend(:, :, :, i_voxel);
    score.Kend(i_voxel) = mean(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    score_std.Kend(i_voxel) = std(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    
    all_scores_nvoxel = all_scores.Riem(:, :, :, i_voxel);
    score.Riem(i_voxel) = mean(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    score_std.Riem(i_voxel) = std(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    
    all_scores_nvoxel = all_scores.Frob(:, :, :, i_voxel);
    score.Frob(i_voxel) = mean(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    score_std.Frob(i_voxel) = std(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    
    all_scores_nvoxel = all_scores.CKA(:, :, :, i_voxel);
    score.CKA(i_voxel) = mean(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
    score_std.CKA(i_voxel) = std(all_scores_nvoxel(~isnan(all_scores_nvoxel)));
end

%% ===== Plotting =====
Fig = figure();
hold on

plot(n_voxels, score.Riem, 'b', 'LineWidth', 2)
plot(n_voxels, score.Frob, 'k', 'LineWidth', 2)
plot(n_voxels, score.Spea, 'r', 'LineWidth', 2)
plot(n_voxels, score.Pear, 'g', 'LineWidth', 2)
plot(n_voxels, score.Kend, 'm', 'LineWidth', 2)
plot(n_voxels, score.CKA, 'c', 'LineWidth', 2)

n_score = n_run*n_ROI*(n_ROI - 1)/2 - 1;

plt = fill([n_voxels, fliplr(n_voxels)],[score.Riem + score_std.Riem/sqrt(n_score), fliplr(score.Riem - score_std.Riem/sqrt(n_score))], 'b', 'EdgeColor', 'none');
alpha(plt, 0.2)
plt = fill([n_voxels, fliplr(n_voxels)],[score.Frob + score_std.Frob/sqrt(n_score), fliplr(score.Frob - score_std.Frob/sqrt(n_score))], 'k', 'EdgeColor', 'none');
alpha(plt, 0.2)
plt = fill([n_voxels, fliplr(n_voxels)],[score.Spea + score_std.Spea/sqrt(n_score), fliplr(score.Spea - score_std.Spea/sqrt(n_score))], 'r', 'EdgeColor', 'none');
alpha(plt, 0.2)
plt = fill([n_voxels, fliplr(n_voxels)],[score.Pear + score_std.Pear/sqrt(n_score), fliplr(score.Pear - score_std.Pear/sqrt(n_score))], 'g', 'EdgeColor', 'none');
alpha(plt, 0.2)
plt = fill([n_voxels, fliplr(n_voxels)],[score.Kend + score_std.Kend/sqrt(n_score), fliplr(score.Kend - score_std.Kend/sqrt(n_score))], 'm', 'EdgeColor', 'none');
alpha(plt, 0.2)
plt = fill([n_voxels, fliplr(n_voxels)],[score.CKA + score_std.CKA/sqrt(n_score), fliplr(score.CKA - score_std.CKA/sqrt(n_score))], 'c', 'EdgeColor', 'none');
alpha(plt, 0.2)

% 2^{nd}-moment
title('2^{nd}-moment', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('# Response channels', 'FontSize', 12)
ylabel('score', 'FontSize', 12)
legend('Riemannian', 'Frobenius', 'Spearman', 'Pearson', 'Kendall''s \tau', 'CKA', 'FontSize', 12, 'Location','southeast')
% legend('Riemannian', 'Frobenius', 'Spearman', 'Pearson', 'Kendall''s \tau', 'FontSize', 12, 'Location','southeast')


ylim([0.5, 1])
xlim([10, inf])
% set(gcf,'Position',[100 100 300 250])

% ===== Save Plot =====
% saveas(gcf, sprintf('results/response_channels_limitation_%s.svg', method))
print(Fig,sprintf('results/response_channels_limitation_%s', method),'-dpdf','-r1000')


%%
function R = get_representation(B, method)
    switch method
        case 'RSM'
            R = corr(B', 'type', 'Pearson');
        case 'G'
            R = secondMoment(B, false);
        case 'G-centered'
            R = secondMoment(B, true);
    end
end