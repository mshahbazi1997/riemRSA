function plot_geodesic(A, B, cl)
    AB = RiemannLogMap(A, B);

    n_t = 1e3;
    t = linspace(0, 1, n_t);
    ab_t = nan(3, n_t);
    mask = [true, false; true, true];
    for i_t = 1 : n_t
        AB_i = RiemannExpMap(A, AB*t(i_t));
        ab_t(:, i_t) = AB_i(mask);
    end
    
    plot3(ab_t(1, :), ab_t(2, :), ab_t(3, :), 'Color', [0.5, 0.5, 1], 'LineWidth', 4, 'Color', cl)
    
    xlabel('(1, 1)', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('(2, 1)', 'FontSize', 20, 'FontWeight', 'bold')
    zlabel('(2, 2)', 'FontSize', 20, 'FontWeight', 'bold')
    
    d = distance_riemann(A, B);
    n_half = floor(n_t/2);
    text(ab_t(1, n_half) + 0.2, ab_t(2, n_half) + 0, ab_t(3, n_half)+ 0.1, sprintf('%.2f', d), 'FontSize', 18, 'Color', cl)

    view([-45, 45])
end