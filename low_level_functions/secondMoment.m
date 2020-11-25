function G = secondMoment(beta, isCKA)
% beta = randn(10,100);
% beta is [#condition x #voxels]
% test ==> centeredG(randn(10,100))
n_voxels = size(beta, 2);
if isCKA
    beta = beta - mean(beta);
end
G = beta*beta';
G = G/n_voxels;
end