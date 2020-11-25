function plot_colorful(x, y, z)
    hold on
    n_t = length(x);
    
    for i_t = 1 : n_t
        if nargin > 2
            plot3(x(i_t), y(i_t), z(i_t), 'x', 'Color', [0.8, 0.8, 0.8]*i_t/n_t)
        else
            plot(x(i_t), y(i_t), 'x', 'Color', [0.8, 0.8, 0.8]*i_t/n_t)
        end
    end
end