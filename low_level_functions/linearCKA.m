function CKA = linearCKA(G_1, G_2)

% G_i is a [#condition x #condition] 
% beta_1 = randn(10,100);
% beta_2 = randn(10,100);
% 
% G_1 = secondMoment(beta_1,1);
% G_2 = secondMoment(beta_2,1);
% for test ===> linearCKA(G_1, G_2)

vecG_1 = reshape(G_1,1,[]);
vecG_2 = reshape(G_2,[],1);

CKA = (vecG_1*vecG_2)/sqrt((vecG_1*vecG_1')*(vecG_2'*vecG_2));
end
