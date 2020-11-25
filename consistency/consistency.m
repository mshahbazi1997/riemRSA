clear
clc

%% ===== Settings =====
subs = [1, 2, 3, 4, 5];

use_rsm = false;

low_level_func_path = ['..', filesep, 'low_level_functions'];
utils_path = ['..', filesep, 'utils'];

dataset_name = 'Kamitani';
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
n_comb_of_runs = size(comb_runs, 1);
do_perm = true;
n_perm = 1e2;

%% ==== Big Loop ====
bias_corrected = nan(data_spec.n_roi, ...
                     n_comb_of_runs, ...
                     data_spec.n_sub, ...
                     n_rep_similarity);
            
h = waitbar(0, 'Please wait...');

for i_sub = 1 : data_spec.n_sub
    fprintf('===> Proccessing subject %d ... \n', subs(i_sub));
    
    for i_roi = 1 : data_spec.n_roi
        waitbar( (i_roi + (i_sub - 1)*data_spec.n_roi)/data_spec.n_sub/data_spec.n_roi, h)
        
        for i_comb = 1 : n_comb_of_runs
            r_1 = comb_runs(i_comb, 1);
            r_2 = comb_runs(i_comb, 2);
            
            b_r_1 = u_sub_run_roi{i_sub, r_1, i_roi};
            b_r_2 = u_sub_run_roi{i_sub, r_2, i_roi};
            
            if use_rsm
                rg_r_1 = corr(b_r_1', 'type', 'Pearson');  % rep. geometry
                rg_r_2 = corr(b_r_2', 'type', 'Pearson');
            else
                rg_r_1 = b_r_1*b_r_1';
                rg_r_2 = b_r_2*b_r_2';
            end
            
            [~, ~, ~, ~, bias_corrected(i_roi, i_comb, i_sub, 1)] = pval_perm_Spearman(rg_r_1, rg_r_2, do_perm, n_perm, use_rsm);
            [~, ~, ~, ~, bias_corrected(i_roi, i_comb, i_sub, 2)] = pval_perm_Pearson(rg_r_1, rg_r_2, do_perm, n_perm, use_rsm);
            [~, ~, ~, ~, bias_corrected(i_roi, i_comb, i_sub, 3)] = pval_perm_Riemann(rg_r_1, rg_r_2, do_perm, n_perm, use_rsm);
            [~, ~, ~, ~, bias_corrected(i_roi, i_comb, i_sub, 4)] = pval_perm_Frobenius(rg_r_1, rg_r_2, do_perm, n_perm, use_rsm);
            [~, ~, ~, ~, bias_corrected(i_roi, i_comb, i_sub, 5)] = pval_perm_Kendall(rg_r_1, rg_r_2, do_perm, n_perm, use_rsm);
            if ~use_rsm
                [~, ~, ~, ~, bias_corrected(i_roi, i_comb, i_sub, 6)] = pval_perm_CKA(rg_r_1, rg_r_2, do_perm, n_perm, use_rsm);
            end
        end
    end
end

close(h)

%% ==== Save Results ====
save([save_path, filesep, 'bias_corrected_', int2str(use_rsm), '.mat'], 'bias_corrected')

