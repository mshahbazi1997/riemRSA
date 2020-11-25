function [u_sub_run_roi, spec] = load_data(name, dataset_folder_path, subs)
    switch name
        case 'Haxby'
            % Set specs
            if nargin < 3
                subs = 1 : 5;
            end
            spec.n_sub = length(subs);
            
            spec.n_run = 11;
            
            important_rois = 81:90;
            spec.n_roi = length(important_rois);
            
            spec.n_cond = 8;
            
            % Load roi labels
            rois_load_addr =  [dataset_folder_path, filesep, 'Haxby2001', filesep, ...
                               'juelich_label_all.nii.gz'];
            roi_labels = niftiread(rois_load_addr);
            
            u_sub_run_roi = cell(spec.n_sub, spec.n_run, spec.n_roi);
            
            for i_sub = 1 : spec.n_sub
                fprintf('loading subject %d data ... \n', i_sub)
                for i_run = 1 : spec.n_run
                    for i_roi = 1 : spec.n_roi
                        roi = important_rois(i_roi);
                        u_sub_run_roi{i_sub, i_run, i_roi} = nan(spec.n_cond, sum(roi_labels(:) == roi));
                    end
                    
                    for i_cond = 1 : spec.n_cond
                        u_row = niftiread([dataset_folder_path, filesep, 'Haxby2001', filesep, ...
                                           's', int2str(i_sub), filesep, ...
                                           'run', int2str(i_run), filesep, ...
                                           'zstat', int2str(i_cond), '_stan.nii.gz']);
                                   
                        for i_roi = 1 : spec.n_roi
                            roi = important_rois(i_roi);
                            u_sub_run_roi{i_sub, i_run, i_roi}(i_cond, :) = ...
                                u_row(roi_labels == roi);
                        end
                    end
                end
            end
            
        case 'Kamitani'
            labels = load([dataset_folder_path, filesep, 'Kamitani', filesep, 'labels']).labels;
            roi_betas = load([dataset_folder_path, filesep, 'Kamitani', filesep, 'ROIbetas']).data_testPercept;
            
            if nargin < 3
                subs = 1 : 5;
            end
            spec.n_sub = length(subs);
            spec.n_run = 35;
            spec.n_roi = 10;
            spec.n_cond = 50;
            
            u_sub_run_roi = cell(spec.n_sub, spec.n_run, spec.n_roi);
            
            for i_sub = 1 : spec.n_sub
                sub = subs(i_sub);
                
                idx_cond_run = nan(spec.n_cond, spec.n_run);
                for i_cond = 1 : spec.n_cond
                    idx_cond_run(i_cond, :) = find(labels{sub} == i_cond);
                end
                
                for i_run = 1 : spec.n_run
                    for i_roi = 1 : spec.n_roi
                        u_sub_run_roi{i_sub, i_run, i_roi} = ...
                            roi_betas{i_sub}{i_roi}(idx_cond_run(:, i_run), :);
                    end
                end
            end
    end
end