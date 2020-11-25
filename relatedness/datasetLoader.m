function all_copes = datasetLoader(sub, Data_base_addr, isBeta)

% constants
n_runs = 11;
n_conds = 8;

% check if we should load beta pattern or z-stat pattern
if isBeta
    fileprefix = 'cope';
else
    fileprefix= 'zstat';
end

all_copes = cell(n_runs, 1);

returnHear = pwd;
cd(['..' filesep '..'])

for run = 1 : n_runs
    copes_load_addr = [Data_base_addr filesep 'sub'...
        sprintf('%.2d', sub) filesep 'run' sprintf('%.2d', run) '.feat' filesep 'stats'];

    cope = niftiread([copes_load_addr, filesep, sprintf([fileprefix '%d_stan.nii'], 1)]);
    copes = nan(numel(cope), n_conds);
    copes(:, 1) = cope(:);

    for cope_i = 2 : n_conds
        cope = niftiread([copes_load_addr, filesep, sprintf([fileprefix '%d_stan.nii'], cope_i)]);
        copes(:, cope_i) = cope(:);
    end
    all_copes{run} = copes;

    fprintf('subject%02d run%02d loaded successfully!\n', sub, run);
end
cd(returnHear)
end