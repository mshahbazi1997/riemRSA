clc; clear;

low_level_func_path = ['..', filesep, 'low_level_functions'];
addpath(low_level_func_path);

addpath('MultivarCon-selected')
% constants
n_layer = 9;
n_trial = 10;
n_sample = 20;
n_metrics = 4;
n_giving_input = 10;


%% calculating score and similarity matrices

comb_trials = combnk(1:n_trial,2);
perm_trials = [comb_trials;comb_trials(:,1),comb_trials(:,2)];

score_maps = cell(n_giving_input,1);


h = waitbar(0, 'Please wait...');

for i_giving_input = 1 : n_giving_input
    
    waitbar(i_giving_input/n_giving_input,h)
        
    % data of layers
    load(['d_list_' num2str(i_giving_input) '.mat'])
    d_list{1} = l1;d_list{2} = l2;d_list{3} = l3;d_list{4} = l4;d_list{5} = l5;
    d_list{6} = l6;d_list{7} = l7;d_list{8} = l8;d_list{9} = l9;
    clear l1 l2 l3 l4 l5 l6 l7 l8 l9

    score_average_across_perm_trials = zeros(n_layer,n_layer,n_metrics);

    for i_pt = 1 : length(perm_trials)

        trial_1 = perm_trials(i_pt, 1);
        trial_2 = perm_trials(i_pt, 2);

        score_average_across_samples = zeros(n_layer,n_layer,n_metrics);

        for i_sample = 1 : n_sample
            similarity = zeros(n_layer,n_layer,n_metrics);
            for lr_1 = 1 : n_layer
                for lr_2 = 1 : n_layer
                    X{1} = d_list{lr_1}(:, :, i_sample, trial_1);
                    Y{1} = d_list{lr_2}(:, :, i_sample, trial_2);

                    % 1: RCA % 2: Riem 3: CKA 4: dCor
                    [similarity(lr_1, lr_2, 1), ~, similarity(lr_1, lr_2, 2), similarity(lr_1, lr_2, 3)] = data2rc(X,Y,'Correlation');
                    similarity(lr_1, lr_2, 4) = data2dCor(X,Y);
                end
            end
            score = zeros(n_layer,n_layer,n_metrics);

            % for RCA
            [~, argmax] = max(similarity(:,:,1), [], 2);
            for lr = 1 : n_layer;score(lr, argmax(lr),1)=1; end
            % for Riem
            [~, argmin] = min(similarity(:,:,2), [], 2);
            for lr = 1 : n_layer;score(lr, argmin(lr),2)=1; end
            % for CKA
            [~, argmax] = max(similarity(:,:,3), [], 2);
            for lr = 1 : n_layer;score(lr, argmax(lr),3)=1; end
            % for dCor
            [~, argmax] = max(similarity(:,:,4), [], 2);
            for lr = 1 : n_layer;score(lr, argmax(lr),4)=1; end

            score_average_across_samples = score_average_across_samples + ...
                score/n_sample;
        end
        score_average_across_perm_trials = score_average_across_perm_trials + ...
            score_average_across_samples/length(perm_trials);
    end
    score_maps{i_giving_input} = score_average_across_perm_trials;
    clear d_list
end
save('score_maps.mat','score_maps')

close(h)

% fprintf('RCA  accuracy: %.3f\n',mean(diag(score_average_across_perm_trials(:,:,1))))
% fprintf('Riem accuracy: %.3f\n',mean(diag(score_average_across_perm_trials(:,:,2))))
% fprintf('CKA  accuracy: %.3f\n',mean(diag(score_average_across_perm_trials(:,:,3))))
% fprintf('dCor accuracy: %.3f\n',mean(diag(score_average_across_perm_trials(:,:,4))))

%% === load ===
load('score_maps.mat')
scores = [];
for i_giving_input = 1 : n_giving_input
    score = [];
    for i_metric = 1 : n_metrics 
        score = cat(1, score, mean(diag(score_maps{i_giving_input}(:,:,i_metric))));
    end
    scores = cat(2, scores, score);
end

%% confusion matrix
scores_confusion_map=zeros(n_layer,n_layer,n_metrics);
for i_giving_input=1:n_giving_input
    for i_metric=1:n_metrics 
        scores_confusion_map(:,:,i_metric)=scores_confusion_map(:,:,i_metric)+(score_maps{i_giving_input}(:,:,i_metric))/n_giving_input;
    end
end

FIG=figure('Color','w');


labels = {'RCA', 'dRiem', 'CKA', 'dCor'};
titles = {'A.','B.','C.','D.'};
order=[4 1 2 3];

for i_metric=1:n_metrics 
    subplot(2,2,i_metric)
    imagesc(100*scores_confusion_map(:,:,order(i_metric)))
    colorbar
    caxis([0 100])
    colormap pink
    axis square;
    title([titles{i_metric} labels{order(i_metric)}])
    xticks(1:2:9);
    yticks(1:2:9)
    
    if i_metric==1 || i_metric==3 
        ylabel('layer number', 'FontSize', 12)
    end
    if i_metric==3 || i_metric==4
        xlabel('layer number', 'FontSize', 12)
    end
end

F = getframe(FIG);
imwrite(F.cdata,fullfile('Graphics','confusion_mats.png'))


%% Bargraph
Fig = figure;
hold on
spread = std(scores, [], 2);%/sqrt(n_giving_input);
meanvl = mean(scores, 2);
c = categorical(labels);
bar(c,meanvl,0.4,'FaceColor',[0.5 0.5 1],'FaceAlpha',1)
errorbar(c,meanvl,spread,'ko','MarkerSize',1,'CapSize',5,'LineWidth',2)
ylabel('Mean +/- SD')
temp = get(gca,'YLim');set(gca,'YLim',[temp(1),temp(2)+.1])
title('Accuracy of identifying corresponding layers')
print(Fig,fullfile('Graphics',['deep_part1']),'-dpdf','-r1000')







