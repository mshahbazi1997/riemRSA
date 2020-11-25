clear; clc; close all
%% === Constants ===
sub = 1;
isBeta = 1;
n_Runs = 11;
n_ROIs = 10;
n_voxels = 1000;
voxels = 10:2:50;
Data_base_addr = ['datasets' filesep 'Haxby2001'];
ROIs_load_addr =  [Data_base_addr filesep 'juelich_label_all.nii.gz'];
%% === Load dependecies ===

low_level_func_path = ['..', filesep, 'low_level_functions'];
addpath(low_level_func_path);

addpath('rsatoolbox-selected');


load(['dependencies' filesep 'design.mat'], '-ascii')
all_copes = datasetLoader(sub, Data_base_addr, isBeta);
%% Simulation options
simulationOptions.X = design(:, 1:2:16); % extract conditions
simulationOptions.nPatterns = 20;
nPatterns = simulationOptions.nPatterns;
%% === Save results ===
save_addr = 'Result';
if ~exist(save_addr,'file')
    mkdir(save_addr)
end
%% === Simulation ===

for ROI = 1:n_ROIs
    fileName = sprintf('CKA_relatedness_sub0%d_ROI%d.mat',sub,ROI);
    
    acc.RSM.spea = nan(n_Runs, length(voxels));
    acc.RSM.pear = nan(n_Runs, length(voxels));
    acc.RSM.kend = nan(n_Runs, length(voxels));
    acc.RSM.riem = nan(n_Runs, length(voxels));
    acc.RSM.frob = nan(n_Runs, length(voxels));
    acc.G.spea = nan(n_Runs, length(voxels));
    acc.G.pear = nan(n_Runs, length(voxels));
    acc.G.kend = nan(n_Runs, length(voxels));
    acc.G.riem = nan(n_Runs, length(voxels));
    acc.G.frob = nan(n_Runs, length(voxels));
    acc.G.CKA = nan(n_Runs, length(voxels));

    h = waitbar(0, sprintf('relatedness sub0%d ROI%d.mat',sub,ROI));
    for run = 1 : n_Runs
        rng((ROI-1)*run+run)
        waitbar(run/n_Runs, h)
        simulationOptions.B_true = generateBetaPatterns(covPattern(all_copes{run}, ROI, ROIs_load_addr), n_voxels);
        simulationOptions.cov = covPattern(all_copes{run}, ROI, ROIs_load_addr);
        % p-values
        p.RSM.spea = nan(nPatterns, length(voxels));
        p.RSM.pear = nan(nPatterns, length(voxels));
        p.RSM.kend = nan(nPatterns, length(voxels));
        p.RSM.riem = nan(nPatterns, length(voxels));
        p.RSM.frob  = nan(nPatterns, length(voxels));
        p.G.spea = nan(nPatterns, length(voxels));
        p.G.pear = nan(nPatterns, length(voxels));
        p.G.kend = nan(nPatterns, length(voxels));
        p.G.riem = nan(nPatterns, length(voxels));
        p.G.frob  = nan(nPatterns, length(voxels));
        p.G.CKA = nan(nPatterns, length(voxels));

        for l = 1 : length(voxels)
            simulationOptions.voxel = voxels(l);
            % OUTPUT is a (nNoisyPatterns x 2 x 11) matrix
            OUTPUT = simulatefMRIData(simulationOptions);

            p.RSM.spea(:, l) = OUTPUT(:, 1, 1);
            p.RSM.pear(:, l) = OUTPUT(:, 1, 2);
            p.RSM.kend(:, l) = OUTPUT(:, 1, 3);
            p.RSM.riem(:, l) = OUTPUT(:, 1, 4);
            p.RSM.frob(:, l) = OUTPUT(:, 1, 5);
            p.G.spea(:, l) = OUTPUT(:, 1, 6);
            p.G.pear(:, l) = OUTPUT(:, 1, 7);
            p.G.kend(:, l) = OUTPUT(:, 1, 8);
            p.G.riem(:, l) = OUTPUT(:, 1, 9);
            p.G.frob(:, l) = OUTPUT(:, 1, 10);
            p.G.CKA(:, l) = OUTPUT(:, 1, 11);
        end

        for l = 1 : length(voxels)

            acc.RSM.spea(run, l) = FDRthreshold(p.RSM.spea(:, l), 0.01, 1)/nPatterns;
            acc.RSM.pear(run, l) = FDRthreshold(p.RSM.pear(:, l), 0.01, 1)/nPatterns;
            acc.RSM.kend(run, l) = FDRthreshold(p.RSM.kend(:, l), 0.01, 1)/nPatterns;
            acc.RSM.riem(run, l) = FDRthreshold(p.RSM.riem(:, l), 0.01, 1)/nPatterns;
            acc.RSM.frob(run, l)  = FDRthreshold(p.RSM.frob(:, l), 0.01, 1)/nPatterns;
            acc.G.spea(run, l) = FDRthreshold(p.G.spea(:, l), 0.01, 1)/nPatterns;
            acc.G.pear(run, l) = FDRthreshold(p.G.pear(:, l), 0.01, 1)/nPatterns;
            acc.G.kend(run, l) = FDRthreshold(p.G.kend(:, l), 0.01, 1)/nPatterns;
            acc.G.riem(run, l) = FDRthreshold(p.G.riem(:, l), 0.01, 1)/nPatterns;
            acc.G.frob(run, l)  = FDRthreshold(p.G.frob(:, l), 0.01, 1)/nPatterns;
            acc.G.CKA(run, l) = FDRthreshold(p.G.CKA(:, l), 0.01, 1)/nPatterns;
        end
    end
    close(h)
    % === Save results ===
    save([save_addr, filesep, fileName], 'acc')
end
