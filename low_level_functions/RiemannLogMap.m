function Out   = RiemannLogMap(P, X)
if (min(eig(X)) <= 0 )
    warning('min eig is %e, X RSM shrinked!', min(eig(X)))
    X = Riemann_shrink(X);
end

[U, Delta] = eig(P);
G = U*sqrt(Delta);
Y = G\X*inv(G)';
[V, Sigma] = eig(Y);
Out = (G*V)*diag(log(diag(Sigma)))*(G*V)';
