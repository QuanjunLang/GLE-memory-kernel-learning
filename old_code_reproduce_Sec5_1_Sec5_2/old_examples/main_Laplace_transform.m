clc
close all
clear all

folder = fileparts(which(mfilename));
addpath(genpath(folder));
%% generate data
N = 5001;
L = 25;
dx = L/(N-1);
xgrid = linspace(dx, L, N);
nsr = 0;

h_grid = mlf(3/2, 1, -abs(xgrid).^(3/2), 10);
h_grid = h_grid + randn(size(xgrid))*nsr;


g_grid = mlf(3/2, 0, -abs(xgrid).^(3/2), 10)./xgrid;
g_grid(1) = (h_grid(2) - h_grid(1))/dx;
g_grid_diff = gradient(h_grid, dx);


gamma_kernel = @(x) (1/sqrt(pi) * x.^(-1/2)).*(x>0);
gamma_grid = gamma_kernel(xgrid);

%% using change of variable
G = zeros(N, N);
for i = 1:N
    for j = 1:N
        G(i,j) = 1/(xgrid(i) + xgrid(j));
    end
end

%% Laplace transform of data
Lap_mat = zeros(N, N);
for i = 1:N
    for j = 1:N
        Lap_mat(i, j) = exp(-xgrid(i)*xgrid(j));
    end
end
h_lap = Lap_mat*h_grid'*dx;
g_lap = Lap_mat*g_grid'*dx;

gh_lap = -Lap_mat*(g_lap./h_lap)*dx;


%%
figure;hold on;
plot(log10(xgrid), G*gamma_grid'*dx);
plot(log10(xgrid), gh_lap);
legend('True b', 'Lap b')

%%

% c = (G + 0.000001*eye(N))\gh_lap
[lambda_opt, c] = L_curve(G, gh_lap, 'RKHS', 0);

figure;hold on;
plot(xgrid, c)
plot(xgrid, gamma_grid)