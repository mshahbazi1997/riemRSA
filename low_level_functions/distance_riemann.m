function a = distance_riemann(A,B)
% Low-level function (not documented yet, see pval_perm_Riemann)
	a = sqrt(sum(log(eig(A,B)).^2));