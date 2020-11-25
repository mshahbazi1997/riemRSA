function OUTPUT = simulatefMRIData(simulationOptions)

    % ==== Simulation ====
    nPatterns = simulationOptions.nPatterns;
    B_true = simulationOptions.B_true;

    voxel = simulationOptions.voxel;
    OUTPUT = nan(nPatterns, 2, 11);
    cov = simulationOptions.cov;

    for o = 1:nPatterns
        
        rng(100+o)
        B_noisy = generateBetaPatterns(cov, voxel);
        output = permutationMeasures(B_true, B_noisy);
        OUTPUT(o, :, :) = output;

    end
end