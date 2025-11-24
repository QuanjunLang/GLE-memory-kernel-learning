clear all
close all
clc

%%
lb = 0;
rb = 10;
knot_num = 10;
free_bdry = 0;
deg = 3;

if free_bdry
    knots = [lb*ones(1, deg-1), linspace(lb, rb, knot_num+1), rb*ones(1, deg-1)];
    N = knot_num+deg-1;
else
    knots = linspace(lb, rb, knot_num+1);
    N = knot_num-deg+1;
end



basis = cell(N, 1);
for i = 1:N
    coef = zeros(N, 1)';
    coef(i) = 1;
    sp = spmak(knots, coef);
    basis{i} = @(x) fnval(sp, x);
end



basis_deg_1 = cell(N-1, 1);
for i = 1:N
    coef = zeros(N, 1)';
    coef(i) = 1;
    sp = spmak(knots, coef);
    basis{i} = @(x) fnval(sp, x);
end



basis_derivative = cell(N, 1);
for i = 1:N-1
    coef_1 = zeros(N, 1);
    coef_1(i) = deg/(knots(i+deg-1) - knots(i));
    
    coef_2 = zeros(N, 1);
    coef_2(i+1) = deg/(knots(i+deg) - knots(i+1));
    
    sp_1 = spmak(knots, coef_1);
    sp_2 = spmak(knots, coef_2);
    
    basis_derivative{i} = @(x) fnval(sp_1, x) - fnval(sp_2, x);
end



i = 5;

xgrid = linspace(lb, rb, 1000);
dx = xgrid(2) - xgrid(1);
f = basis{i}(xgrid);
df = gradient(f, dx);
g = basis_derivative{i}(xgrid);

figure;hold on;
plot(xgrid, f);
plot(xgrid, df);
plot(xgrid, g);


%%
tau = linspace(0,pi,101); k = 4;
knots = augknt(linspace(0,pi,11),k);
colmat = spcol(knots,k,brk2knt(tau,3));
colmat = spcol(knots,k,brk2knt(tau,3));
plot(colmat(2:3:end,:))
% coefs = (colmat(3:3:end,:)/5-colmat(1:3:end,:))\(-sin(2*tau).');
% sp = spmak(knots,coefs.');


%%
lb = 0;
rb = 10;
knot_num = 20;
deg = 3;
basis_derivative = 3;
xgrid = linspace(0, 10, 101);
free_bdry = 1;
[basis_mat, knots] = spline_basis_new(lb, rb, knot_num, deg, basis_derivative, xgrid, free_bdry)

