clc
close all
clear all
%% generate data
N = 101;
L = 10;
dx = L/(N-1);
xgrid = linspace(dx, L, N);

% h = @(x) exp(-x).*sin(x*2);
h_grid = mlf(3/2, 1, -abs(xgrid).^(3/2), 10);
% g should be the derivative of h
g_grid = mlf(3/2, 0, -abs(xgrid).^(3/2), 10)./xgrid;
g_grid(1) = (h_grid(2) - h_grid(1))/dx;
% h_grid = sin(xgrid);
% h = @(x) sin(x);

% plot(g_grid);hold on;plot(gradient(h_grid)/dx);
%% kernel K using function handle h
% K = zeros(N, N);
% for i = 1:N
%     for j = 1:N
%         x = xgrid(i);
%         y = xgrid(j);
%         K(i, j) = sum(h(xgrid - x).*h(xgrid - y).*(xgrid >x).*(xgrid > y))*dx;
%     end
% end
% K1 = K;
%% Using direct integral
% h_pad = [h_grid, zeros(size(xgrid))];
% K = zeros(N, N);
% t_ind = 1:N;
% for i = 1:N
%     for j = 1:N
%         k = max(i, j);
%         t_ind = k:N;
%         t_i_ind = t_ind - i + 1;
%         t_j_ind = t_ind - j + 1;
%         K(i, j) = sum(h_grid(t_i_ind).*h_grid(t_j_ind))*dx;
%     end
% end
% K2 = K;
%% using change of variable
h_pad = [h_grid, zeros(size(xgrid))];
KK = zeros(N+1, 1);
h_pad = [h_grid, zeros(size(xgrid))];
for t = 1:N+1
    KK(t) = sum(h_grid.*h_pad(t:t+N-1))*dx;
end

K = zeros(N, N);
for i = 1:N
    for j = 1:N
        K(i, j) = KK(abs(i-j)+1);
    end
end

%% generate cross term
g_pad = [g_grid, zeros(size(xgrid))];
LL = zeros(N, 1);
for i = 1:N
    LL(i) = -sum(h_grid.*g_pad(i:N+i-1))*dx;
end


%%
gamma_kernel = @(x) (1/sqrt(pi) * x.^(-1/2)).*(x>0);
gamma_est = K\LL;
plot(xgrid, gamma_kernel(xgrid));hold on;
plot(xgrid, 20*gamma_est);

%% generate g using direct integration
g_direct = zeros(N, 1);
for i = 1:N
    t = xgrid(i);
    sgrid = xgrid (1:i);
    g_direct(i) = sum(gamma_kernel(t - sgrid).*h_grid(1:i))*dx;
end


%%
% figure;
% subplot(121)
% plot(xgrid, h_grid)
% subplot(122)
% 
% [XX, YY] = meshgrid(xgrid, xgrid);
% surf(XX, YY, K);shading flat;

%%
% figure;
% plot(xgrid, )