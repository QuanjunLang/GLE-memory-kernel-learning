clc
close all
clear all

%%
N = 1000;

f = @(x) exp(-1*x) + 1./(1+10*x).^2;
xgrid = 10.^linspace(-10, 10, N);
% xgrid = linspace(0.00001, 10, N);

plot(log10(xgrid), f(xgrid))
ylim([0, 3])

%%
G = @(x, y) 1./(x+y);

A = zeros(N, N);
for i = 1:N
    for j = 1:N
        A(i, j) = G(xgrid(i), xgrid(j));
    end
end


%%
[XX, YY] = meshgrid(log10(xgrid), log10(xgrid));
surf(XX, YY, A)
shading flat

%%
[U, S, V] = svd(A);

%%
figure;plot(U(:, 1:5))
figure;plot(log10(diag(S)))



%%
sgrid = linspace(-10, 10, 1000);
plot(sgrid, pi./cosh(pi*sgrid))