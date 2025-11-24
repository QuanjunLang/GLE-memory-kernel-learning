clc
close all
clear all
rng(1)

%%
a = 1;
R = @(x) exp(-abs(x)*a);
S = @(w) 1/pi*a./(w.^2 + a^2);
wu = 10*pi;
N = 800;
M = 2^13;
f = generate_correlated_noise(R, S, wu, N, M)


%%
R = @(x) (1-3*x.^2)./((1 + x.^2).^3);
S = @(w) 1/4*w.^2.*exp(-abs(w));

wu = 4*pi;
N = 128;
M = 2560;
dt = 0.2;
f = generate_correlated_noise(R, S, wu, N, M)