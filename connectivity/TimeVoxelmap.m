close all
clear;
clc;

% choose the nonlinearity type
% nonlinearity = 1 ===> absolute value function 
% nonlinearity = 2 ===> relu function
% nonlinearity = 3 ===> square function
% nonlinearity = 4 ===> linear function


% Parameter settings
nSubj = 50;         % number of subjects (replications)
nRuns = 2;          % number of runs. MVPD requires at least 2 runs.
sigma = 1;          % std of Noise in ROI1 and ROI2

low_level_func_path = ['..', filesep, 'low_level_functions'];
addpath(low_level_func_path);

addpath('MultivarCon-selected')


% Optional parameters
opt.zscore = 0; % Whether to Z-score timeseries for dCor (need to turn off for Examples 5-6)
opt.nRandomisation = 100; % Needs to be >1 to get more stable estimates of null distributions (though will trigger parfor if >1)

max_nVox = 150;
min_nVox = 50;
n_samples = 10;
nVox = floor(linspace(min_nVox, max_nVox, n_samples));


for nonlinearity = 4:4

    if nonlinearity == 1
        filename = 'abs';
    elseif nonlinearity == 2
        filename = 'relu';
    elseif nonlinearity == 3
        filename = 'square';
    else
        filename = 'linear';
    end

    % Simulation

    for v = 1 : length(nVox)
        nTimes = floor(linspace(0.1*nVox(v), 0.7*nVox(v), n_samples));
        for t = 1 : length(nTimes)

            nVoxs = [nVox(v) nVox(v)];    % number of voxels in ROI1 and ROI2
            nTime = nTimes(t);                % number of time points
            T = randn(nVoxs); 
            C = eye(nVoxs(1)); % correlation within ROI
            X = {}; Y = {};

            for s=1:nSubj
                for r=1:nRuns
                    X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);

                    if nonlinearity == 1
                        Y{s}{r} = abs(X{s}{r}*T);
                    elseif nonlinearity == 2
                        X{s}{r} = max(X{s}{r},0);
                        Y{s}{r} = max(X{s}{r}*T,0);
                    elseif nonlinearity == 3
                        Y{s}{r} = abs(X{s}{r}*T).^2;
                    else
                        Y{s}{r} = X{s}{r}*T;
                    end
                    X{s}{r} = X{s}{r} + sigma*randn(nTime,nVoxs(1)); % Add independent measurement noise
                    Y{s}{r} = Y{s}{r} + sigma*randn(nTime,nVoxs(2));
                end
            end
            [MVconn(v,t),MVconn_null(v,t)] = computeMVconn(X,Y,opt);
        end
        clc
        disp(v)
    end
    save(['matFile', filesep, filename '.mat'], 'MVconn', 'MVconn_null')
    clear MVconn MVconn_null
end
