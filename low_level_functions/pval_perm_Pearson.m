function [p, pe_corr, perm_corrs, z, bias_corrected_norm] = pval_perm_Pearson(A, B, do_perm, varargin)
    % [p, dist, perm_dists, z, bias_corrected_dist] = pval_perm_Pearson(A, B, do_perm, n_rept, isRSM)
    %
    % Calculates the Pearson correlation of two representational matrices
    % inputs:
    %   - A: representational matrix: [#cond x #cond]
    %   - B: representational matrix: [#cond x #cond], permuted one.
    %   - do_perm: if true, condition permutation will be done.
    %   - n_rept (optional, default is 1000): number of permutations to be done
    %   - isRSM (optional, default is true)
    % outputs:
    %   - pe_corr: Pearson correlation of A and B
    %  if dp_perm is true, following will be returned. o.w. they will be
    %  nan:
    %   - p: p-value of permutation test.
    %   - perm_corrs: all permuted correlations: [n_rept x 1]
    %   - z: t-stat*sqrt(n_perm)
    %   - bias_corrected_dist
    
    %% ======== handling inputs ========
    p = inputParser;
    addRequired(p, 'A');
    addRequired(p, 'B');
    addRequired(p, 'do_perm');
    addOptional(p, 'n_rept', 1e3);
    addOptional(p, 'isRSM', true);
    
    parse(p, A, B, do_perm, varargin{:}) 
    
    A = p.Results.A;
    B = p.Results.B;
    do_perm = p.Results.do_perm;
    n_rept = p.Results.n_rept;
    isRSM = p.Results.isRSM;
    %% =================================
    
    vectorize = @(X) X(triu(ones(size(X,1)), isRSM)==1);
    L = size(A, 1);
    
    pe_corr = corr(vectorize(A), vectorize(B), 'type', 'Pearson');    
    
    % permutation
    if do_perm
        perm_corrs = nan(n_rept, 1);
        for perm_i = 1 : n_rept
            rnd_perm = randperm(L);
            perm_corrs(perm_i) = corr(vectorize(A), vectorize(B(rnd_perm, rnd_perm)), 'type', 'Pearson');
        end
        
        % tstat & p-value calculation
        p = 0; % ranksum(perm_corrs, pe_corr, 'tail', 'left', 'method', 'exact'); 
        
        mu = mean(perm_corrs);
        sigma = std(perm_corrs);
        z = (pe_corr - mu)/sigma;
        bias_corrected_norm = pe_corr - mu;
    else
        p = 0;
        z = 0;
        perm_corrs = nan;
        bias_corrected_norm = 0;
    end

end