function Out   = RiemannExpMap(P,X)

[U, Delta] = eig(P);
G = U*sqrt(Delta);
Y = G\X*inv(G)';
[V, Sigma] = eig(Y);
Out = (G*V)*diag(exp(diag(Sigma)))*(G*V)';