function [MVconn,MVconn_null] = computeMVconn(X,Y,opt)
% this is mainly a wrapper function: given the input MV patterns for 2
% regions, it computes MVPD, LPRD, dCor, RC and MIM and also FC and FC_SVD
% the function also simulates the null hypothesis of no functional
% connectivity by randomly shuffling the time points (first component) in
% each run of every subject. This assumes no temporal autocorrelation in
% data and would not be appropriate for temporally smooth data.
% results from permuted data are returned in MVconn_null, which contains
% one null number per subject per iteration. The number of iterations are
% specified in opt.nRandomisation.
%
% inputs:
%        X: a cell array, one entry per subject. X{s} is the data from
%        subject "s" and contains mutiple cell arrays, one per run. For
%        example X{2}{1} is the data from run1 of subject2. The data is 
%        time x voxels.
%        Y : the same as X for the second region
% 
% outputs:
%        MVconn: contains the following fields:
%               FC, FCSVD, FCCCA, MVPD, LPRD, dCor, RCA (MIM, ImCohSVD, MVLagCoh, LagCohSVD
%               each is a column vector with one number per subject.
%        MVconn_null: contains the same fields as MVconn. Each would be a
%        matrix with size of nSubjects x nRandomisations.
% Hamed Nili


if ~isfield(opt,'nRandomisation')
    opt.nRandomisation = 1;
end
if ~isfield(opt,'zscore')
    opt.zscore = 0;
end

nSub = length(X);

% zscore patterns within each run if opt.zscore is set in the options
if opt.zscore
    for r = 1:numel(X)
        X{r} = zscore(X{r},0,2);
        Y{r} = zscore(Y{r},0,2);
    end
end

% Calculate connectivity on given data
for s=1:nSub
    [dcor(s,1),~] = data2dCor(X{s},Y{s});
    [rc(s,1),rc_riem_rsm(s,1), rc_riem_G(s,1), rc_cka(s,1)] = data2rc(X{s},Y{s},'Correlation');
end

% Calculate connectivity when X and Y independent random noise (since
% some connectivity measures, eg dCor, not bounded by 0 or -1)

bdcor       = NaN(nSub,opt.nRandomisation);
brc         = NaN(nSub,opt.nRandomisation);
brc_riem_G      = NaN(nSub,opt.nRandomisation);
brc_riem_rsm    = NaN(nSub,opt.nRandomisation);
brc_cka         = NaN(nSub,opt.nRandomisation);

if opt.nRandomisation == 1 %if only one, then don't bother with parfor like below
    if nSub < 20
        warning('May not be sufficient subjects/randomisations to estimate null properly')
    end
    iter = 1;
    for s=1:nSub % Ensure reasonably accurate estimate
        bX = {}; bY = {};
        for r=1:length(X{s})
            bX{r} = X{s}{r}(randperm(size(X{s}{r},1)),:);
            bY{r} = Y{s}{r}(randperm(size(Y{s}{r},1)),:);
        end
        [bdcor(s,iter),~] = data2dCor(bX,bY);
        [brc(s,iter),brc_riem_rsm(s,iter), brc_riem_G(s,iter),brc_cka(s,iter)] = data2rc(bX,bY,'Correlation');
    end
elseif opt.nRandomisation > 1
    for s=1:nSub % Ensure reasonably accurate estimate
        fprintf('null subject %d from %d \n',s,nSub)
        parfor iter = 1:opt.nRandomisation
            bX = {}; bY = {};
            for r=1:length(X{s})
                bX{r} = X{s}{r}(randperm(size(X{s}{r},1)),:);
                bY{r} = Y{s}{r}(randperm(size(Y{s}{r},1)),:);
            end
            [tmp_bdcor{iter},~] = data2dCor(bX,bY);
            [tmp_brc{iter},tmp_brc_riem_rsm{iter}, tmp_brc_riem_G{iter},tmp_brc_cka{iter}] = data2rc(bX,bY,'Correlation');
        end

        bdcor(s,:) = cat(2,tmp_bdcor{:});
        brc(s,:) = cat(2,tmp_brc{:});
        brc_riem_G(s,:) = cat(2,tmp_brc_riem_G{:});
        brc_riem_rsm(s,:) = cat(2,tmp_brc_riem_rsm{:});
        brc_cka(s,:) = cat(2,tmp_brc_cka{:});
    end
end
fprintf('\n')


MVconn.dCor = dcor;
MVconn.RCA = rc;
MVconn.RCARiemG = rc_riem_G;
MVconn.RCARiemRSM = rc_riem_rsm;
MVconn.RCACKA = rc_cka;

MVconn_null.dCor = mean(bdcor,2);
MVconn_null.RCA = mean(brc,2);
MVconn_null.RCARiemG = mean(brc_riem_G,2);
MVconn_null.RCARiemRSM = mean(brc_riem_rsm,2);
MVconn_null.RCACKA = mean(brc_cka,2);

return

