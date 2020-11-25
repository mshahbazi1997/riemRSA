clear
clc
returnHere = pwd; % We'll come back here later

%% ===== User input =====
n_cond = 8;
n_sub = 20;

%% ===== Add paths =====
depend_path = 'dependencies';

toolbox_path = 'rsatoolbox-selected';
addpath(genpath(toolbox_path));

low_level_funcs_path = '../low_level_functions';
addpath(low_level_funcs_path)

%% ===== Options =====
userOptions.voxelSize = [3, 3, 3.75];
userOptions.searchlightRadius = 15;
userOptions.RDMcorrelationType = 'all';
userOptions.analysisName = '12th_Aban';
userOptions.rootPath = 'simulation_junks';
userOptions.ModelColor = [0, 1, 0];

searchlightOptions.averageSessions = true;
searchlightOptions.monitor = false;
searchlightOptions.nConditions = n_cond;
searchlightOptions.nSessions = 1;
searchlightOptions.fisher = false;

simulationOptions = simulationOptions_demo_SL();

%% ===== Load model =====
models = constructModelRDMs(modelRDMs_SL_sim, userOptions);

%% ===== Load dependencies =====
load([depend_path, filesep, 'sampleMask_org.mat'])
load([depend_path, filesep, 'anatomy.mat'])
         
%% ===== Simulation =====
for subI = 1 : n_sub
    fprintf('simulating fullBrain volumes for subject %d \n', subI)

    subject = ['subject', num2str(subI)];
    maskName = 'mask';
    
    [B_true, Mask, Y_true, fMRI_sub] = simulateClusteredfMRIData_fullBrain(simulationOptions);
    B_noisy = fMRI_sub.B;
    singleSubjectVols = B_noisy';
    
    mask = m;
    
    fprintf('computing correlation and distance maps for subject %d \n',subI)
    [smm_dists, smm_corrs, ns] = ...
        searchlightMapping_fMRI(singleSubjectVols, models, mask, userOptions, searchlightOptions);
    
    gotoDir('results', 'maps');
    save(['result_',subject,'.mat'], 'smm_dists', 'smm_corrs');
    
    clear smm_dists smm_corrs searchlightRDMs
    cd(returnHere)
end
