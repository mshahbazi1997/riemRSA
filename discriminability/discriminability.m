clear
clc

%% ===== Settings =====
subs = [1, 2, 3, 4, 5];

use_rsm = true;

low_level_func_path = ['..', filesep, 'low_level_functions'];
utils_path = ['..', filesep, 'utils'];

dataset_name = 'Haxby';
data_load_path = ['..', filesep, '..', filesep, 'datasets'];

save_path = ['results', filesep, dataset_name];
mkdir(save_path);

%% ===== Add paths =====
addpath(low_level_func_path);
addpath(utils_path);

%% ===== Load data =====
[u_sub_run_roi, data_spec] = load_data(dataset_name, data_load_path, subs);

fprintf('Data loaded successfully!\n');

%% ===== Constants =====
if use_rsm
    n_rep_similarity = 5;
else
    n_rep_similarity = 6;
end
comb_runs = combnk(1 : data_spec.n_run, 2);
comb_rois = combnk(1 : data_spec.n_roi, 2);
n_comb_of_runs = size(comb_runs, 1);
n_comb_of_rois = size(comb_rois, 1);
do_perm = true;
n_perm = 1e2;

%% ===== Inline function =====
within_minus_between = @(func, rg_run_1_roi_1, rg_run_2_roi_1, rg_run_1_roi_2, rg_run_2_roi_2) ...
    out5(func, rg_run_1_roi_1, rg_run_2_roi_1, do_perm, n_perm, use_rsm) ...
    + out5(func, rg_run_1_roi_2, rg_run_2_roi_2, do_perm, n_perm, use_rsm) ...
    - out5(func, rg_run_1_roi_1, rg_run_2_roi_2, do_perm, n_perm, use_rsm) ...
    - out5(func, rg_run_2_roi_1, rg_run_1_roi_2, do_perm, n_perm, use_rsm);

%% ===== Big Loop =====
bias_corrected = nan(data_spec.n_roi, ...
                     data_spec.n_roi, ...
                     n_comb_of_runs, ...
                     data_spec.n_sub, ...
                     n_rep_similarity);
            
h = waitbar(0, 'Please wait...');

for i_sub = 1 : data_spec.n_sub
    fprintf('===> Proccessing subject %d ... \n', subs(i_sub));
    
    for i_comb_roi = 1 : n_comb_of_rois
        waitbar( (i_comb_roi + (i_sub - 1)*n_comb_of_rois)/data_spec.n_sub/n_comb_of_rois, h)
        
        roi_1 = comb_rois(i_comb_roi, 1);
        roi_2 = comb_rois(i_comb_roi, 2);
        
        if roi_1 == roi_2
            continue;
        end
        
        for i_comb_run = 1 : n_comb_of_runs
            run_1 = comb_runs(i_comb_run, 1);
            run_2 = comb_runs(i_comb_run, 2);
            
            if run_1 == run_2
                continue;
            end
            
            b_run_1_roi_1 = u_sub_run_roi{i_sub, run_1, roi_1};
            b_run_2_roi_1 = u_sub_run_roi{i_sub, run_2, roi_1};
            b_run_1_roi_2 = u_sub_run_roi{i_sub, run_1, roi_2};
            b_run_2_roi_2 = u_sub_run_roi{i_sub, run_2, roi_2};
            
            if use_rsm
                rg_run_1_roi_1 = corr(b_run_1_roi_1', 'type', 'Pearson');  % rep. geometry
                rg_run_2_roi_1 = corr(b_run_2_roi_1', 'type', 'Pearson');
                rg_run_1_roi_2 = corr(b_run_1_roi_2', 'type', 'Pearson');
                rg_run_2_roi_2 = corr(b_run_2_roi_2', 'type', 'Pearson');
            else
                rg_run_1_roi_1 = b_run_1_roi_1*b_run_1_roi_1';
                rg_run_2_roi_1 = b_run_2_roi_1*b_run_2_roi_1';
                rg_run_1_roi_2 = b_run_1_roi_2*b_run_1_roi_2';
                rg_run_2_roi_2 = b_run_2_roi_2*b_run_2_roi_2';
            end
            
            bias_corrected(roi_1, roi_2, i_comb_run, i_sub, 1) = ...
                within_minus_between(@pval_perm_Spearman, rg_run_1_roi_1, rg_run_2_roi_1, rg_run_1_roi_2, rg_run_2_roi_2);
            bias_corrected(roi_1, roi_2, i_comb_run, i_sub, 2) = ...
                within_minus_between(@pval_perm_Pearson, rg_run_1_roi_1, rg_run_2_roi_1, rg_run_1_roi_2, rg_run_2_roi_2);
            bias_corrected(roi_1, roi_2, i_comb_run, i_sub, 3) = ...
                within_minus_between(@pval_perm_Riemann, rg_run_1_roi_1, rg_run_2_roi_1, rg_run_1_roi_2, rg_run_2_roi_2);
            bias_corrected(roi_1, roi_2, i_comb_run, i_sub, 4) = ...
                within_minus_between(@pval_perm_Frobenius, rg_run_1_roi_1, rg_run_2_roi_1, rg_run_1_roi_2, rg_run_2_roi_2);
            bias_corrected(roi_1, roi_2, i_comb_run, i_sub, 5) = ...
                within_minus_between(@pval_perm_Kendall, rg_run_1_roi_1, rg_run_2_roi_1, rg_run_1_roi_2, rg_run_2_roi_2);
            
            if ~use_rsm
                bias_corrected(roi_1, roi_2, i_comb_run, i_sub, 6) = ...
                    within_minus_between(@pval_perm_CKA, rg_run_1_roi_1, rg_run_2_roi_1, rg_run_1_roi_2, rg_run_2_roi_2);
            end
        end
    end
end

close(h)

%% ==== Save Results ====
save([save_path, filesep, 'bias_corrected_', int2str(use_rsm), '.mat'], 'bias_corrected')

