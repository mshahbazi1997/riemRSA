function [p, fro_norm, perm_dists, z, bias_corrected_norm] = pval_perm_Frobenius(A, B, do_perm, varargin)
    % [p, dist, perm_dists, z, bias_corrected_dist] = pval_perm_Frobenius(A, B, do_perm, n_rept, isRSM)
    %
    % Calculates the Frobenius norm of difference of two representational matrices
    % inputs:
    %   - A: representational matrix: [#cond x #cond]
    %   - B: representational matrix: [#cond x #cond], permuted one.
    %   - do_perm: if true, condition permutation will be done.
    %   - n_rept (optional, default is 1000): number of permutations to be done
    %   - isRSM: useless (for compatability)
    % outputs:
    %   - dist: Frobenius norm of difference of A and B
    %  if dp_perm is true, following will be returned. o.w. they will be
    %  nan:
    %   - p: p-value of permutation test.
    %   - perm_dists: all permuted distances: [n_rept x 1]
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
    %% =================================
    
    L = size(A, 1);
    
    fro_norm = norm(A - B, 'fro');    
    % permutation
    if do_perm
        perm_dists = nan(n_rept, 1);
        for perm_i = 1 : n_rept
            rnd_perm = randperm(L);
            perm_dists(perm_i) = norm(A - B(rnd_perm, rnd_perm), 'fro');
        end
        
        % tstat & p-value calculation
        p = 0; % ranksum(perm_dists, fro_norm, 'tail', 'right', 'method', 'exact'); 
        
        mu = mean(perm_dists);
        sigma = std(perm_dists);
        z = (mu - fro_norm)/sigma;
        bias_corrected_norm = mu - fro_norm;
    else
        p = 0;
        z = 0;
        bias_corrected_norm = 0;
        perm_dists = nan;
    end

end