close all
clear all
clc

%% 
dx = 0.02;
L = 200;
xgrid = (0:dx:L)';

h = @(x) (x+1).^(-2);


d = @(x) exp(-1*x)/2;

h_grid = h(xgrid);
d_grid = d(xgrid(2:end));
H = tril(toeplitz((h_grid(1:end-1) + h_grid(2:end))/2));


Hd = H*d_grid*dx;

figure;hold on;
plot(xgrid(2:end), d_grid);
plot(xgrid, h_grid);
plot(xgrid(2:end), Hd);

error = sum(Hd.^2)*dx
d_norm = sum(d_grid.^2)*dx
h_norm = sum(h_grid.^2)*dx

d_norm*h_norm