if plot_1d_ON
    all_v_sorted = sortrows(all_v')';
    all_dv_sorted = gradient(all_v_sorted', dt)';
    figure;
    bd = max(abs(all_v_sorted), [], 'all');
    [aa, bb]=size(all_v);
    K = min(500, bb);
    c = linspace(1,10,K);
    for i = 1:500
        % scatter(all_v_sorted(i, 1:K), all_dv_sorted(i, 1:K), [], c, 'filled')
        scatter(all_v_sorted(i, 1:K), all_dv_sorted(i, 1:K)*0, [], c, 'filled')
        ylim([-0.5, 0.5])
        xlim([-bd, bd])
        pause(0.001)
    end


end


if plot_2d_ON
    all_v_sorted = sortrows(all_v')';
    all_dv_sorted = gradient(all_v_sorted', dt)';
    figure;
    bd = max(abs([all_v_sorted, all_dv_sorted]), [], 'all');

    K = 500;
    c = linspace(1,10,K);
    for i = 1:500
        scatter(all_v_sorted(i, 1:K), all_dv_sorted(i, 1:K), [], c, 'filled')
        % scatter(all_v_sorted(i, 1:K), all_dv_sorted(i, 1:K)*0+1, [], c, 'filled')
        xlim([-bd, bd])
        ylim([-bd, bd])
        pause(0.0000001)
    end
end

