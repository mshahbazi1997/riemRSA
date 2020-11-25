function RSM_shr = Riemann_shrink(RSM)
% Low-level function (not documented yet, see pval_perm_Riemann)

    max_num_iter = 3;
    shr_ratio = 0.1;
    
    RSM_shr = RSM;
    while max_num_iter > 0
        RSM_shr = RSM_shr + shr_ratio*eye(size(RSM));
        if min(eig(RSM_shr)) > 0
            break;
        end
    end
    
end