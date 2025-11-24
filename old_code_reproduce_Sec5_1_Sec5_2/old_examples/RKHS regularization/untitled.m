close all
clear all
clc

%%
n = 2;
syms x
w = sym('w', [n, 1]);
z = sym('z', [n, 1]);

h = 0;
for i = 1:n
    h = h + w(i)/(x - z(i));
end

g = 1/h;
%%
pretty(g)
simp = partfrac(g, 'FactorMode', 'complex');
% ss = subs(simp, w, ones(n, 1)/2)
pretty(simp)

%%
C = children(simp)