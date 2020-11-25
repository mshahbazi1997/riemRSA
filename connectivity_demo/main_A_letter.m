clear
clc
close all

%%
low_level_func_path = ['..', filesep, 'low_level_functions'];
addpath(low_level_func_path);

%%

n_p = 50;
sig = 0.1;
t = linspace(0, 1, n_p);
s_y = [2*t, -2*t, 2*t - 1] + sig*randn(1, 3*n_p);
s_x = -[-2*(2*t) + 2, 2*(-2*t) + 2, 0*t] + sig*randn(1, 3*n_p);
s_z = zeros(size(s_x));

%%
R = roty(-45);
xyz = [s_x', s_y', s_z']*R;

r = 3;

x = xyz(:, 1) + r;
y = xyz(:, 2);
z = xyz(:, 3) + r;

figure
plot_colorful(x, y, z);

xlabel('(1, 1)', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('(1, 2)', 'FontSize', 18, 'FontWeight', 'bold')
zlabel('(2, 2)', 'FontSize', 18, 'FontWeight', 'bold')

xlim([0, 2*r])
zlim([0, 2*r])
ylim([-r, r])

view([45, 45])

grid on

%%
figure

subplot(2, 2, 1)
plot_mds(x, y, z, 'Riem');

subplot(2, 2, 2)
plot_mds(x, y, z, 'Frob');

%
%A = [1, 0; ...
%     -0.6, 1];
A = rand(2);
x_p = nan(size(x));
y_p = nan(size(y));
z_p = nan(size(z));

for i_t = 1 : length(x)

    R_t = [x(i_t), y(i_t); ...
           y(i_t), z(i_t)];

    R_t_p = A*R_t*A';

    x_p(i_t) = R_t_p(1, 1);
    y_p(i_t) = R_t_p(2, 1);
    z_p(i_t) = R_t_p(2, 2);
end

%
subplot(2, 2, 3)
plot_mds(x_p, y_p, z_p, 'Riem');

subplot(2, 2, 4)
plot_mds(x_p, y_p, z_p, 'Frob');

%%
figure
imagesc(A)

axis square

caxis([0, 1])

set(gca,'visible','off')
set(gca,'Xtick',[])