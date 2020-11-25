function uv = plot_mds(x, y, z, meter)
    n_t = length(x);

    D = zeros(n_t, n_t);

    for i_t_1 = 1 : n_t
        for i_t_2 = i_t_1 + 1 : n_t
            
            R_1 = [x(i_t_1), y(i_t_1); ...
                   y(i_t_1), z(i_t_1)];
            R_2 = [x(i_t_2), y(i_t_2); ...
                   y(i_t_2), z(i_t_2)];
            
            switch meter
                case 'Frob'
                    D(i_t_1, i_t_2) = norm(R_1 - R_2, 'fro');
                case 'Riem'
                    D(i_t_1, i_t_2) = real(distance_riemann(R_1, R_2));
            end
        end
    end
    
    %% MDS
    D = D + D';
    
    uv = mdscale(D, 2);

    plot_colorful(uv(:, 1), -uv(:, 2))
    
    set(gca,'visible','off')
    set(gca,'Xtick',[])
end