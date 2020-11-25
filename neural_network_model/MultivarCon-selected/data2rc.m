function [rc_w,rc_b, rc_riem_G, rc_cka] = data2rc(X,Y,metric,corrType)
% This function computes the representational connectivity from the
% multivariate time series of two regions. 
% It computes both the wihtin- and between-run RDM correlations for the two
% regions (rc_w and rc_b respectively).
% Hamed Nili



nVoxs = zeros(1, 2);
[nTimes, nVoxs(1)] = size(X{1});
[~ , nVoxs(2)] = size(Y{1});


if nargin < 3
    metric = 'correlation';
end
if ~exist('corrType','var')
    corrType = 'Pearson';
end

nruns = numel(X); % number of runs

for r=1:nruns
    rdms1(:,r)    = pdist(X{r},metric)';
    rdms2(:,r)    = pdist(Y{r},metric)';
end
corrmat  = corr(rdms1,rdms2,'type',corrType);
rc_w = mean(diag(corrmat));
rc_b = mean(corrmat(logical(1-eye(nruns))));



%% ====  Riemannian ====
if (nTimes/min(nVoxs))<=1
    
    for r=1:nruns
        Gs1(:, :, r) = secondMoment(X{r}, 0);
        Gs2(:, :, r) = secondMoment(Y{r}, 0);
    end

%     for i = 1:nruns
%         
%         if (min(eig(Gs1(:, :, r))) <=eps )
%             Gs1(:, :, r) = Riemann_shrink(Gs1(:, :, r));
%             warning('Gs1 shrinked!')
%         end
% 
%         if (min(eig(Gs2(:, :, r))) <=eps )
%             Gs2(:, :, r) = Riemann_shrink(Gs2(:, :, r));
%             warning('rdm2 matrix shrinked!')
%         
%         end
%     end

    dist_riem_G  = zeros(nruns,1);
    for i = 1:nruns
        dist_riem_G(i) = distance_riemann(Gs1(:, :, i), Gs2(:, :, i));
    end
    rc_riem_G = mean(dist_riem_G);

else
    
    mVoxs = min(nVoxs);
    batch_size = floor(mVoxs/2);
    
    
    for i_batch = 1 : floor(nTimes/batch_size) + 1
        if i_batch <= floor(nTimes/batch_size)
            for r=1:nruns
                dist_riem_G(r, i_batch) = distance_riemann(secondMoment(X{r}((i_batch-1)*batch_size+1:(i_batch)*batch_size, :), 0),...
                                                                                       secondMoment(Y{r}((i_batch-1)*batch_size+1:(i_batch)*batch_size, :), 0));
            end
        else
            if nTimes>(i_batch-1)*batch_size
                for r=1:nruns
                    dist_riem_G(r, i_batch) = distance_riemann(secondMoment(X{r}((i_batch-1)*batch_size+1:end, :), 0),...
                                                                                           secondMoment(Y{r}((i_batch-1)*batch_size+1:end, :), 0));
                end
            end
        end
    end
    rc_riem_G = mean(mean(dist_riem_G));
%     rc_riem_G = nan;
end
% 
% 
%% ==== CKA ====

clear Gs1
clear Gs2

for r=1:nruns
    Gs1(:, :, r) = secondMoment(X{r}, 1);
    Gs2(:, :, r) = secondMoment(Y{r}, 1);
end

% riem_distances = zeros(nruns, 1);
if nruns > 1
    ckas  = zeros(nruns,1);
    for i = 1:nruns
        ckas(i) = linearCKA(Gs1(:, :, i), Gs2(:, :, i));
    end
    rc_cka = mean(ckas);
else
    rc_cka = linearCKA(Gs1(:, :, 1), Gs2(:, :, 1));
end


