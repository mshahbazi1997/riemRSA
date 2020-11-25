clear; clc
addpath('low level functions')

% This script is only for ploting the resut

nonlinearity = 2;

if nonlinearity == 1
    filename = 'abs';
elseif nonlinearity == 2
    filename = 'relu';
elseif nonlinearity == 3
    filename = 'square';
else 
    filename = 'linear';
end

FIG=figure('name',filename,'Color','w','Position',[1 1 370 2*330]);


load(['matFile', filesep, filename '.mat'])


% constants
max_nVox = 150;
min_nVox = 50;
n_samples = 10;
nVox = floor(linspace(min_nVox, max_nVox, n_samples));

methods = {'dCor','RCA','RCARiemG', 'RCACKA'};
titles = {'A. ','B. ','C. ','D. ','E. '};


%% color map 2
comb_methods =  combnk(1:numel(methods), 2);%combnk(1:numel(methods), 2);

map = zeros(n_samples, n_samples, 3);
txt = cell([size(map, 1), size(map, 2)]);
labels_short = {'D', 'r', 'R', 'C'};

for v=1:n_samples
    for t = 1:n_samples
        scores = zeros(numel(methods), numel(methods));
        for c = 1 : length(comb_methods)
            if (strcmp(methods{comb_methods(c,1)},'RCARiemG') || strcmp(methods{comb_methods(c,2)},'RCARiemG') )
                if strcmp(methods{comb_methods(c,1)},'RCARiemG') 
                    allval = - MVconn(v,t).(methods{comb_methods(c,1)});
                    allnull = - MVconn_null(v,t).(methods{comb_methods(c,1)});
                    left = (allval-allnull)/std(allval-allnull);
                    
                    allval = MVconn(v,t).(methods{comb_methods(c,2)});
                    allnull = MVconn_null(v,t).(methods{comb_methods(c,2)});
                    right = (allval-allnull)/std(allval-allnull);
                    
                    [pval, ~, ~] = signrank(left, right, 'tail', 'right', 'method', 'approximate');
                    
                    if pval < 0.05
                        scores(comb_methods(c,1), comb_methods(c,2)) = 1;
                        scores(comb_methods(c,2), comb_methods(c,1)) = -1;
                    elseif pval > 0.95
                        scores(comb_methods(c,1), comb_methods(c,2)) = -1;
                        scores(comb_methods(c,2), comb_methods(c,1)) = 1;
                    else
                        scores(comb_methods(c,1), comb_methods(c,2)) = 0;
                        scores(comb_methods(c,2), comb_methods(c,1)) = 0;
                    end
                    
                else
                    allval = MVconn(v,t).(methods{comb_methods(c,1)});
                    allnull = MVconn_null(v,t).(methods{comb_methods(c,1)});
                    left = (allval-allnull)/std(allval-allnull);
                    
                    allval = -MVconn(v,t).(methods{comb_methods(c,2)});
                    allnull = -MVconn_null(v,t).(methods{comb_methods(c,2)});
                    right = (allval-allnull)/std(allval-allnull);
                    
                    [pval, ~, ~] = signrank(left, right, 'tail', 'right', 'method', 'approximate');
                    
                    if pval < 0.05
                        scores(comb_methods(c,1), comb_methods(c,2)) = 1;
                        scores(comb_methods(c,2), comb_methods(c,1)) = -1;
                    elseif pval > 0.95
                        scores(comb_methods(c,1), comb_methods(c,2)) = -1;
                        scores(comb_methods(c,2), comb_methods(c,1)) = 1;
                    else
                        scores(comb_methods(c,1), comb_methods(c,2)) = 0;
                        scores(comb_methods(c,2), comb_methods(c,1)) = 0;
                    end
                end
            else
                allval = MVconn(v,t).(methods{comb_methods(c,1)});
                allnull = MVconn_null(v,t).(methods{comb_methods(c,1)});
                left = (allval-allnull)/std(allval-allnull);

                allval = MVconn(v,t).(methods{comb_methods(c,2)});
                allnull = MVconn_null(v,t).(methods{comb_methods(c,2)});
                right = (allval-allnull)/std(allval-allnull);

                [pval, ~, ~] = signrank(left, right, 'tail', 'right', 'method', 'approximate');

                if pval < 0.05
                    scores(comb_methods(c,1), comb_methods(c,2)) = 1;
                    scores(comb_methods(c,2), comb_methods(c,1)) = -1;
                elseif pval > 0.95
                    scores(comb_methods(c,1), comb_methods(c,2)) = -1;
                    scores(comb_methods(c,2), comb_methods(c,1)) = 1;
                else
                    scores(comb_methods(c,1), comb_methods(c,2)) = 0;
                    scores(comb_methods(c,2), comb_methods(c,1)) = 0;
                end
            end    
        end
        
        
        scores = sum(scores,2);
        ind_max = find(scores == max(scores));
        if length(ind_max) == 1
            map(v,t, :) = index2color(ind_max);
            txt{v, t} = labels_short{ind_max};
        else
            col = zeros(1,3);
            txt{v, t} = '';
            for idx = 1 : length(ind_max)
                col = col + index2color(ind_max(idx));
%                 map(v,t, ind_max(idx)) = 1/length(ind_max);
%                 txt{v, t} = [txt{v, t}, ',', labels_short{ind_max(idx)}];
            end
            map(v,t, :) = col;%1/length(ind_max);
            
%             txt{v, t} = txt{v, t}(2:end);
        end
    end
end


subplot(2, 1, 2)
hold on
imshow(map, 'InitialMagnification', 'fit')
axis image;
axis on;

mt = title([titles{end} 'Best-performing measure map']);

set(mt,'Position',[5.5 -0.5],'VerticalAlignment','top','Color',[0 0 0])

xlh = xlabel('$\frac{\#conditions}{\#voxels}$', 'Interpreter', 'latex', 'FontSize', 16);
xlh.Position(2) = xlh.Position(2);  % move the label 0.1 data-units further down
xlh.Position(1) = xlh.Position(1);  % move the label 0.1 data-units further down

ylabel('$\#voxels$', 'Interpreter', 'latex', 'FontSize', 16)
    
cmap = [index2color(1);
        index2color(3);
        index2color(4)];
          
labels = {'dCor', 'RCA', 'dRiem', 'CKA'};
selected_meth = [1 3 4];

% imlegend(cmap ,labels(selected_meth))

xticks([1, 10]);
xticklabels([0.1, 0.7])
yticks([1, 10])
yticklabels([50, 150])

%% ploting each metric's significance map
for v=1:n_samples
    for t = 1:n_samples
        for meth = 1 : length(methods)
            if strcmp(methods{meth},'RCARiemG')
                allval = - MVconn(v,t).(methods{meth});
                allnull = - MVconn_null(v,t).(methods{meth});
                z_maps(v, t, meth) = mean(allval-allnull)/std(allval-allnull);
            else
                allval = MVconn(v,t).(methods{meth});
                allnull = MVconn_null(v,t).(methods{meth});               
                z_maps(v, t, meth) = mean(allval-allnull)/std(allval-allnull);
            end
        end
    end
end

min_min=min(reshape(z_maps(:,:,1),1,[]));
min_max=max(reshape(z_maps(:,:,1),1,[]));
for meth=2:length(methods)
    if min(reshape(z_maps(:,:,meth),1,[]))<min_min
        min_min=min(reshape(z_maps(:,:,meth),1,[]));
    end
    if max(reshape(z_maps(:,:,meth),1,[]))<min_max
        min_max=max(reshape(z_maps(:,:,meth),1,[]));
    end
end

for meth = 1 : length(methods)
%     subplot(3,2,meth)
    subplot(4, 2, meth)

    imagesc(z_maps(:,:,meth))
    axis square;

    colormap gray
%     caxis([1.96 inf])
%     caxis([0 1.645])
%     caxis([0 inf])
    caxis([min_min min_max])
    colorbar
    title([titles{meth} labels{meth}])
    
    xticks([1, 10]);
    xticklabels([0.1, 0.7])
    yticks([1, 10])
    yticklabels([50, 150])
    
    if meth == 3 || meth == 4
        xlabel('$\frac{\#conditions}{\#voxels}$', 'Interpreter', 'latex', 'FontSize', 12)
    end
    
    if meth == 1 || meth == 3
        ylabel('$\#voxels$', 'Interpreter', 'latex', 'FontSize', 12)
    end
end

% saveas(gcf,fullfile('Graphics',[filename '.svg']),'svg')
print(FIG,fullfile('Graphics',filename),'-dpdf','-r1000')


%% function
function c = index2color(index)
    if index == 1
        c = [0 1 0]/1.2;
    elseif index == 2
        c = [0 0 0]/1.2;
    elseif index == 3
        c = [1 0 0]/1.2;
    elseif index == 4
        c= [0 0 1]/1.2;
    end
end

function imlegend(colorArr, labelsArr)
    % For instance if two legend entries are needed:
    % colorArr =
    %   Nx3 array of doubles. Each entry should be an RGB percentage value between 0 and 1
    %
    % labelsArr =
    %   1×N cell array
    %     {'First name here'}    {'Second name here'}    {'etc'}
    hold on;
    for ii = 1:length(labelsArr)
        % Make a new legend entry for each label. 'color' contains a 0->255 RGB triplet
        scatter(10, ii, 1, colorArr(ii,:), 'filled', 'DisplayName', labelsArr{ii});
    end
    hold off;

    legend();
end