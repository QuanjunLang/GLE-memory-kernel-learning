clc
close all
clear all

addpath(genpath(fileparts(which(mfilename))));
%% generate data
I = system_settings();

%% generate K
I = generate_K(I);

%% generate A

I = generate_A(I);
%% generate b
% using integration by parts
I = generate_b(I);

% using integration by parts
I = generate_b_direct(I);

% using direct computation with true gamma
dx = I.dx;
basis = I.basis_mat{1};
basis_derivative = I.basis_mat{2};
H = I.H;
n = I.n;
K = I.K;
gamma_grid = I.gamma_grid;


b = zeros(n, 1);
for i = 1:n
    b(i) = basis(:, i)'*K*gamma_grid'*dx*dx;
end


figure;plot(b);hold on;plot(I.b);plot(I.b_direct);
legend('true gamma b', 'IBP b', 'direct b')
%%
% c = L_curve(I.A, b, 'ID', 1);
% c = (I.A+0.000000001*pinv(I.A))\I.b;
[lambda_opt, c] = L_curve(I.A, I.b_direct, 'RKHS', 1);
sp = spmak(I.knots, c');
f = @(x) fnval(sp, x);

c_noreg = I.A\I.b_direct;
sp = spmak(I.knots, c_noreg');
f_noreg = @(x) fnval(sp, x);
figure;hold on;
plot(I.xgrid, f(I.xgrid))
plot(I.xgrid, I.gamma_grid)
plot(I.xgrid, f_noreg(I.xgrid))

%% generate cross term
dx = I.dx;
basis = I.basis_mat{1};
basis_derivative = I.basis_mat{2};
H = I.H;
n = I.n;
xgrid = I.xgrid;
N = I.N;
h_grid = I.h_grid;
K = I.K;
gamma_grid = I.gamma_grid;


g_grid_diff = gradient(I.h_grid, I.dx);
g_pad = [g_grid_diff, zeros(size(I.xgrid))];
LL = zeros(N, 1);
for i = 1:N
    LL(i) = -sum(h_grid.*g_pad(i:N+i-1))*dx;
end

figure;hold on;
plot(LL);
% plot(LL_1)
plot(K*gamma_grid')
legend('LL', 'true')
%% 
gamma_est = K\LL;
[~, gamma_est_reg] = L_curve(K, LL, 'RKHS', 1);
figure;hold on;
plot(log10(xgrid), gamma_est);
plot(log10(xgrid), gamma_grid);
plot(log10(xgrid), gamma_est_reg);
legend('est', 'true', 'est_reg')

%%
figure;hold on;
plot(xgrid, gamma_est);
plot(xgrid, gamma_grid);
plot(xgrid, gamma_est_reg);
legend('est', 'true', 'est_reg')



% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%% Lap %%%%%%%%%%%%%%%%%
% G = zeros(N, N);
% for i = 1:N
%     for j = 1:N
%         G(i,j) = 1/(xgrid(i) + xgrid(j));
%     end
% end
% 
% %% Laplace transform of data
% Lap_mat = zeros(N, N);
% for i = 1:N
%     for j = 1:N
%         Lap_mat(i, j) = exp(-xgrid(i)*xgrid(j));
%     end
% end
% h_lap = Lap_mat*h_grid'*dx;
% g_lap = Lap_mat*g_grid'*dx;
% 
% gh_lap = -Lap_mat*(g_lap./h_lap)*dx;
% 
% 
% %%
% figure;hold on;
% plot(log10(xgrid), G*gamma_grid'*dx);
% plot(log10(xgrid), gh_lap);
% legend('True b', 'Lap b')
% 
% %%
% 
% % c = (G + 0.000001*eye(N))\gh_lap
% [lambda_opt, c] = L_curve(G, gh_lap, 'RKHS', 1);
% 
% figure;hold on;
% plot(xgrid, c)
% plot(xgrid, gamma_grid)
% 
% 
% 
% %%
% delta = 0;
% HA = K + delta*G;
% Hb = LL + delta*gh_lap;
% 
% % Hc = HA\Hb;
% [lambda_opt, Hc] = L_curve(HA, Hb, 'RKHS', 1);
% 
% figure;hold on;
% plot(xgrid, Hc)
% plot(xgrid, gamma_grid)