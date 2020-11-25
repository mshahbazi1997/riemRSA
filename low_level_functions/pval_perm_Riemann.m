function [p, dist, perm_dists, z, bias_corrected_dist] = pval_perm_Riemann(A, B, do_perm, varargin)
    % [p, dist, perm_dists, z, bias_corrected_dist] = pval_perm_Riemann(A, B, do_perm, n_rept, isRSM)
    %
    % Calculates the Riemannian distance of two representational matrices
    % inputs:
    %   - A: representational matrix: [#cond x #cond]
    %   - B: representational matrix: [#cond x #cond], permuted one.
    %   - do_perm: if true, condition permutation will be done.
    %   - n_rept (optional, default is 1000): number of permutations to be done
    %   - isRSM: useless (for compatability)
    % outputs:
    %   - dist: Riemannian distance of A and B
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
    
    if (min(eig(A)) <=0 )
        A = Riemann_shrink(A);
        warning('A matrix shrinked!')
    end
    
    if (min(eig(B)) <=0 )
        B = Riemann_shrink(B);
        warning('B matrix shrinked!')
    end
    
    dist = distance_riemann(A, B);
    
    % permutation
    if do_perm
        perm_dists = nan(n_rept, 1);
        for perm_i = 1 : n_rept
            rnd_perm = randperm(L);
            perm_dists(perm_i) = distance_riemann(A, B(rnd_perm, rnd_perm));
        end
        
        % tstat & p-value calculation
        p = 0; %ranksum(perm_dists, dist, 'tail', 'right', 'method', 'exact');
        
        mu = mean(perm_dists);
        sigma = std(perm_dists);
        z = (mu - dist)/(sigma);
        bias_corrected_dist = mu - dist;
    else
        p = nan;
        z = nan;
        bias_corrected_dist = nan;
        perm_dists = nan;
    end

end