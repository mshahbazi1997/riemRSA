function [p_map_ttest, bc_maps] = get_p_map_ttest(load_addr, subs, RDMcorrelationType, measure)
    %% Load data
    all_bc_map_1 = load([load_addr, filesep, 'result_subject', int2str(1), '.mat']);
    switch RDMcorrelationType
        case 'corr'
            bc_map_1 = all_bc_map_1.(['smm_', RDMcorrelationType, 's']).([measure, '_bc']);
        case 'dist'
            bc_map_1 = all_bc_map_1.(['smm_', RDMcorrelationType, 's']).([measure, '_bc']);
    end
    
    n_sub = size(subs, 2);
    bc_maps = nan([n_sub, size(bc_map_1)]);
    bc_maps(1, :, :, :) = bc_map_1;
    for i_sub = 1 : n_sub
        sub = subs(i_sub);
        
        all_bc_map_sub = load([load_addr, filesep, 'result_subject', int2str(sub), '.mat']);
        switch RDMcorrelationType
            case 'corr'
                bc_maps(i_sub, :, :, :) = all_bc_map_sub.(['smm_', RDMcorrelationType, 's']).([measure, '_bc']);
            case 'dist'
                bc_maps(i_sub, :, :, :) = all_bc_map_sub.(['smm_', RDMcorrelationType, 's']).([measure, '_bc']);
        end
    end
    
    %% T-test
    [~, p_map_ttest] = ttest(bc_maps, 0, 'Tail', 'right');

    p_map_ttest = squeeze(p_map_ttest);
    
end