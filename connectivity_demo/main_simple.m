clear
clc
close all

%%
low_level_func_path = ['..', filesep, 'low_level_functions'];
addpath(low_level_func_path);

%%
X = wishrnd(eye(2), 3);
Y = wishrnd(eye(2), 3);

A = rand(2);
X_p = A*X*A';
Y_p = A*Y*A';

%%
figure
hold on
view([45, 45])

grid on

plot_mat(X, 'r')
plot_mat(Y, 'r')

plot_mat(X_p, 'r')
plot_mat(Y_p, 'r')

plot_geodesic(X, Y, [0.5, 0.5, 1])
plot_geodesic(X_p, Y_p, [0.5, 0.5, 1])

plot_mat2mat(X, X_p, [0.5, 0.5, 0.5])
plot_mat2mat(Y, Y_p, [0.5, 0.5, 0.5])

