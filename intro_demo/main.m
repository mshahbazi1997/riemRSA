clear
clc
close all

%%
low_level_func_path = ['..', filesep, 'low_level_functions'];
addpath(low_level_func_path);

%%
A = [1, 0.6; 0.6, 1];
B = [1, 0.2; 0.2, 1];
C = [1, 0.9; 0.9, 1];

diff_Riem = distance_riemann(A, B) - distance_riemann(A, C)
diff_Frob = norm(A - B, 'fro') - norm(A - C, 'fro')

%%
AB = RiemannLogMap(A, B);
AC = RiemannLogMap(A, C);

n_t = 1e3;
t = linspace(0, 1, n_t);
ab_t = nan(3, n_t);
ac_t = nan(3, n_t);
mask = [true, false; true, true];
for i_t = 1 : n_t
    AB_i = RiemannExpMap(A, AB*t(i_t));
    AC_i = RiemannExpMap(A, AC*t(i_t));
    ab_t(:, i_t) = AB_i(mask);
    ac_t(:, i_t) = AC_i(mask);
end

%%

Fig=figure('Color','w');
hold on



plot3(ab_t(1, :), ab_t(2, :), ab_t(3, :), 'Color', [0.5, 0.5, 1], 'LineWidth', 3)
plot3(ac_t(1, :), ac_t(2, :), ac_t(3, :), 'Color', [0.5, 0.5, 1], 'LineWidth', 3)

a = A(mask);
b = B(mask);
c = C(mask);
plot3([a(1), b(1)], [a(2), b(2)], [a(3), b(3)], 'Color', [0, 0, 0], 'LineWidth', 2)
plot3([a(1), c(1)], [a(2), c(2)], [a(3), c(3)], 'Color', [0, 0, 0], 'LineWidth', 3)

text(a(1) + 3e-3, a(2), a(3), 'A', 'Color', 'r', 'FontSize', 24) 
text(b(1) + 3e-3, b(2), b(3), 'B', 'Color', 'r', 'FontSize', 24)
text(c(1) + 3e-3, c(2), c(3), 'C', 'Color', 'r', 'FontSize', 24)

plot3(a(1), a(2), a(3), 'ro', 'LineWidth', 20) 
plot3(b(1), b(2), b(3), 'ro', 'LineWidth', 20)
plot3(c(1), c(2), c(3), 'ro', 'LineWidth', 20)

text(ab_t(1, n_t/2), ab_t(2, n_t/2), ab_t(3, n_t/2) - 7e-3, sprintf('%.2f', distance_riemann(A, B)), 'Color', 'b', 'FontSize', 16)
text(ac_t(1, n_t/2), ac_t(2, n_t/2), ac_t(3, n_t/2) - 7e-3, sprintf('%.2f', distance_riemann(A, C)), 'Color', 'b', 'FontSize', 16)

text((a(1) + b(1))/2 + 3e-3, (a(2) + b(2))/2, (a(3) + b(3))/2, sprintf('%.2f', norm(A(mask) - B(mask), 'fro')), 'Color', 'k', 'FontSize', 16)
text((a(1) + c(1))/2 + 3e-3, (a(2) + c(2))/2, (a(3) + c(3))/2, sprintf('%.2f', norm(A(mask) - C(mask), 'fro')), 'Color', 'k', 'FontSize', 16)

axis([0.96, 1.04, 0.05, 0.85, 0.96, 1.04])

xlabel('(1, 1)', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('(2, 1)', 'FontSize', 20, 'FontWeight', 'bold')
zlabel('(2, 2)', 'FontSize', 20, 'FontWeight', 'bold')

view([-45, 45])


% F = getframe(Fig);
% imwrite(F.cdata,fullfile('demo1.pdf'))


print(Fig,'demo1','-dpng','-r1000')
