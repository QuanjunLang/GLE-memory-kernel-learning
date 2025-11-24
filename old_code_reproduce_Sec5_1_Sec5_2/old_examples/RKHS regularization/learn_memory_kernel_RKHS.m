function result = learn_memory_kernel_RKHS(I)

plotON = 1;
%% generate K
I = generate_K(I);

%% generate A

I = generate_A(I);

%% generate b
% using integration by parts
I = generate_b(I);

% using direct computation
I = generate_b_direct(I);

% true b
I = generate_b_true(I);

%% Regression using basis functions
% [~, c_IBP] = L_curve(I.A, I.b, 'RKHS', 0);
% sp = spmak(I.knots, c_IBP');
% gamma_IBP = @(x) fnval(sp, x);
%
% c_IBP_noreg = I.A\I.b;
% sp = spmak(I.knots, c_IBP_noreg');
% gamma_IBP_noreg = @(x) fnval(sp, x);


[~, c_direct] = L_curve(I.A, I.b_direct, 'RKHS', 0);
sp = spmak(I.knots, c_direct');
gamma_direct = @(x) fnval(sp, x);



c_direct_noreg = I.A\I.b_direct;
sp = spmak(I.knots, c_direct_noreg');
gamma_direct_noreg = @(x) fnval(sp, x);


if plotON
    figure;
    subplot(131);hold on;
    plot(I.b_true);
    plot(I.b);
    plot(I.b_direct);
    legend('true gamma b', 'IBP b', 'direct b');
    title('b')
    
    subplot(132);hold on;
    % plot(I.xgrid, gamma_IBP(I.xgrid),           'DisplayName','IBP')
    plot(I.xgrid, gamma_direct(I.xgrid),        'DisplayName','direct')
    % plot(I.xgrid, gamma_IBP_noreg(I.xgrid),     'DisplayName','IBP_noreg')
    plot(I.xgrid, gamma_direct_noreg(I.xgrid),  'DisplayName','direct_noreg')
    plot(I.xgrid, I.gamma(I.xgrid),             'DisplayName','true')
    legend;
    
    subplot(133);hold on;
    % plot(log10(I.log_xgrid), log10(abs(gamma_IBP(I.log_xgrid))),            'DisplayName','IBP')
    plot(log10(I.xgrid), log10(abs(gamma_direct(I.xgrid))),         'DisplayName','direct')
    % plot(log10(I.log_xgrid), log10(abs(gamma_IBP_noreg(I.log_xgrid))),      'DisplayName','IBP_noreg')
    plot(log10(I.xgrid), log10(abs(gamma_direct_noreg(I.xgrid))),   'DisplayName','direct_noreg')
    plot(log10(I.xgrid), log10(I.gamma(I.xgrid)),                   'DisplayName','true')
    legend;
    
    title('Kernel result, using basis')
    
    
end
%% Regression using the grid directly
dx = I.dx;
xgrid = I.xgrid;
N = I.N;
h_grid = I.h_grid;
K = I.K;
gamma_grid = I.gamma_grid;

g_pad = [I.g_grid, zeros(size(I.xgrid))];
LL = zeros(N, 1);
for i = 1:N
    LL(i) = -sum(h_grid.*g_pad(i:N+i-1))*dx;
end



gamma_est = K\LL;
[~, gamma_est_reg] = L_curve(K, LL, 'RKHS', 0);


%%
if plotON
    figure;
    subplot(131);hold on;
    plot(LL);
    plot(K*gamma_grid')
    legend('LL', 'true')
    title('b grid')
    
    subplot(132);hold on;
    plot(xgrid, gamma_est, 'DisplayName','est');
    
    plot(xgrid, gamma_est_reg, 'DisplayName','est_reg');
    plot(xgrid, gamma_grid, 'DisplayName','true');
    legend;
    title('\gamma')
    
    
    subplot(133);hold on;
    plot(log10(xgrid), log10(abs(gamma_est)),'DisplayName', 'est');
    
    plot(log10(xgrid), log10(abs(gamma_est_reg)), 'DisplayName','est_reg');
    plot(log10(xgrid), log10(abs(gamma_grid)), 'DisplayName','true');
    legend;
    title('\gamma, log scale')
end


result.gamma_basis = gamma_direct;
result.gamma_basis_noreg = gamma_direct_noreg;
result.gamma_grid = gamma_est;
result.gamma_grid_noreg = gamma_est_reg;
end
