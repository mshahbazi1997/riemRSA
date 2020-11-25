function Models = modelRDMs_SL_sim()
% Ali's design
% model_mat = kron([0.5, 1; 1, 0.5], ones(4));
% Models.example_1 = model_mat - diag(diag(model_mat));

% Toolbox's default
% Models.main_clusters = kron([
% 			0 0 1 1
% 			0 0 1 1
% 			1 1 0 0
% 			1 1 0 0], ones(10,10));
rdm_struct = load('dependencies/mdl_rdm.mat');
Models.example_1 = rdm_struct.rdm;