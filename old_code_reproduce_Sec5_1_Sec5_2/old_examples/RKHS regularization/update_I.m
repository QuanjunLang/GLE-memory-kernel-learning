function I = update_I(I)


%% grids
I.dx = I.L/(I.N-1);
I.xgrid = linspace(I.dx, I.L, I.N);
[I.basis_mat, I.knots, I.n] = spline_basis_new(I.dx, I.L, I.knot_num, I.deg, I.deg, I.xgrid, I.free_bdry);
I.log_xgrid = 10.^(linspace(-3, 3, I.N));


switch I.example_type
    case 'mlf'
        I.h         = @(x) mlf(3/2, 1, -abs(x).^(3/2), 10);
        I.gamma     = @(x) (1/sqrt(pi) * x.^(-1/2));
        I.Lap_gamma = @(x) x.^(-1/2);
        I.F         = @(x) 0*x;
end

I.gamma_grid = I.gamma(I.xgrid);
switch I.h_data_type
    case 'grid'
        I.h_grid = I.h(I.xgrid) + randn(size(I.xgrid))*I.obs_std;
        
    case 'prony_h'
        I.h_grid = I.prony_h(I.xgrid);
end


I.g_method = 'finite_difference';
switch I.g_method
    case 'analytic'
        g_grid_analytic = mlf(3/2, 0, -abs(I.xgrid).^(3/2), 10)./I.xgrid;
        g_grid_analytic(1) = (I.h_grid(2) - I.h_grid(1))/I.dx;
        I.g_grid = g_grid_analytic;
    case 'finite_difference'
        
        g_grid_diff = gradient(I.h_grid, I.dx);
        I.g_grid = g_grid_diff;
    case 'direct'
        g_grid_direct = zeros(1, I.N);
        for i = 2:I.N
            sgrid = I.xgrid(2:i);
            h_on_s = fliplr(I.h_grid(1:i-1));
            g_grid_direct(i) = -sum(I.gamma(sgrid).*h_on_s)*I.dx;
        end
        I.g_grid = g_grid_direct;
end



I.obs_xgrid = (0:I.obs_dx:(I.prony_N-1)*I.obs_dx)';
I.obs_h_grid = I.h(I.obs_xgrid) + randn(size(I.obs_xgrid))*I.prony_obs_std;








end