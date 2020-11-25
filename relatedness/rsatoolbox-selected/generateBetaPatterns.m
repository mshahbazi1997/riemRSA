function OUT = generateBetaPatterns(Cov, nVoxels)

% Generates a [nConditions nVoxels]-sized beta matrix 
    nConditions = size(Cov, 1);
    mu = zeros(1, nConditions);
    OUT = mvnrnd(mu, Cov,nVoxels)';
end%function