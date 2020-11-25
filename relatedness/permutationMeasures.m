function output = permutationMeasures(B_true, B_noisy)

RSM_true = corr(B_true', 'type', 'Pearson');
RSM_noisy = corr(B_noisy', 'type', 'Pearson');
G_true = secondMoment(B_true, 0);
G_noisy = secondMoment(B_noisy, 0);

centeredG_true = secondMoment(B_true, 1);
centeredG_noisy = secondMoment(B_noisy, 1);

do_perm = true;
n_perm = 1e3;

output = nan(1, 2, 11);

[output(1, 1, 1), ~, ~, ~, output(1, 2, 1)] = pval_perm_Spearman(RSM_true, RSM_noisy, do_perm, n_perm, 1);
[output(1, 1, 2), ~, ~, ~, output(1, 2, 2)] = pval_perm_Pearson(RSM_true, RSM_noisy, do_perm, n_perm, 1);
[output(1, 1, 3), ~, ~, ~, output(1, 2, 3)] = pval_perm_Kendall(RSM_true, RSM_noisy, do_perm, n_perm, 1);
[output(1, 1, 4), ~, ~, ~, output(1, 2, 4)] = pval_perm_Riemann(RSM_true, RSM_noisy, do_perm, n_perm);
[output(1, 1, 5), ~, ~, ~, output(1, 2, 5)] = pval_perm_Frobenius(RSM_true, RSM_noisy, do_perm, n_perm);
[output(1, 1, 6), ~, ~, ~, output(1, 2, 6)] = pval_perm_Spearman(G_true, G_noisy, do_perm, n_perm, 0);
[output(1, 1, 7), ~, ~, ~, output(1, 2, 7)] = pval_perm_Pearson(G_true, G_noisy, do_perm, n_perm, 0);
[output(1, 1, 8), ~, ~, ~, output(1, 2, 8)] = pval_perm_Kendall(G_true, G_noisy, do_perm, n_perm, 0);
[output(1, 1, 9), ~, ~, ~, output(1, 2, 9)] = pval_perm_Riemann(G_true, G_noisy, do_perm, n_perm);
[output(1, 1, 10), ~, ~, ~, output(1, 2, 10)] = pval_perm_Frobenius(G_true, G_noisy, do_perm, n_perm);
[output(1, 1, 11), ~, ~, ~, output(1, 2, 11)] = pval_perm_CKA(centeredG_true, centeredG_noisy, do_perm, n_perm);

end
