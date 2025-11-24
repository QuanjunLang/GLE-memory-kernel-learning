function result_gamma = get_gamma_from_prony_h_esternal_force(I, I_g, result_h, result_g)
%% partial fractional decomposition

w = result_h.w;
lam = result_h.lam;
p = I.prony_p;
h0 = result_h.h(0);
%% Find the inverse of L[h]
syms x
expr = 0;

for i = 1:p
    expr = expr + w(i)/(x - lam(i));
end

expr_inv = vpa(1/expr); 
simp = partfrac(expr_inv, 'FactorMode', 'complex');
% pretty(vpa(simp))

%% Get the laplace transform of memory kernel 

cc = children(simp);
% C = h0*children(simp);

C = zeros(length(cc), 1)+x - x;
for i = 1:length(cc)
    C(i) = cc{i}*h0;
end
% for i = 1:length(C)
%     vpa(C(i))
% end

temp = 0;
w_seq = zeros(length(C) - 2, 1);
lam_seq = zeros(length(C) - 2, 1);

for i = 2:length(C) - 1
    [N, D] = numden(C(i));
    coef = coeffs(D);
    par = coef(2);
    w_seq(i-1) = N/par;
    lam_seq(i-1) = coef(1)/par;
end
lam_seq = -lam_seq;
const = C(end);
%% singularity of gamma ------ negative lambda
for i = 1:length(lam_seq)
    if real(lam_seq(i)) > 0
        lam_seq(i) = -1e-5 + imag(lam_seq(i))*1i;
    end
end

%% extra term f
w_f = result_g.w;
lam_f = result_g.lam;
p_f = length(w_f);

expr_f = 0;

for i = 1:p_f
    expr_f = expr_f + w_f(i)/(x - lam_f(i));
end

expr_f_inv = vpa(1/expr_f); 
% simp_f = partfrac(expr_f_inv, 'FactorMode', 'complex');


a = 1;

%%
% denominator = simp - C(1)/h0;
denominator = simp;
prod_simp = vpa(denominator*expr_f);

% The hard way
% all_terms = expand(prod_simp);
% CC = children(all_terms)
% all_part_frac = 0*x;
% for i = 1:length(CC)
%     i
%     current_term = vpa(CC(i));
%     all_part_frac = all_part_frac + partfrac(CC(i), 'FactorMode', 'complex');
% end
% 
% 
a = 1;


all_part_frac = partfrac(simplify(prod_simp), 'FactorMode', 'complex');

C_f = children(all_part_frac);
for i = 1:length(C_f)
    vpa(C_f(i));
end

w_seq_f = zeros(length(C_f), 1);
lam_seq_f = zeros(length(C_f), 1);

for i = 1:length(C_f)
    [N_f, D_f] = numden(C_f{i});
    coef = coeffs(D_f);
    if length(coef) == 2
        par = coef(2);
        w_seq_f(i) = N_f/par;
        lam_seq_f(i) = coef(1)/par;
    else
        const_f = N_f;
        w_seq_f(i) = [];
        lam_seq_f(i) = [];
    end
    
    a = 1;
end
lam_seq_f = -lam_seq_f;



%%

w_seq = [w_seq;w_seq_f];
lam_seq = [lam_seq;lam_seq_f];

%%
theta_hat = @(x) 0*x + const;
for i = 1:length(w_seq)
    theta_hat = @(x) theta_hat(x) + w_seq(i)./(x - lam_seq(i));
end

%% Inverse Laplace transform (explicit formula)
f = @(x) 0*x;
% f = @(x) const*normpdf(x, 0, 1e-3);
for i = 1:length(w_seq)
    f = @(x) f(x) + w_seq(i)*exp(lam_seq(i).*x);
end




result_gamma.lam_seq = lam_seq;
result_gamma.w_seq = w_seq;
result_gamma.gamma = f;
result_gamma.gamma_hat = theta_hat;


S = @(x) 0*x;
for k = 1:length(w_seq)
    S = @(x) S(x) + w_seq(k)*(-lam_seq(k))/pi ./(x.^2 + lam_seq(k)^2);
end

result_gamma.gamma_fourier = S;


plotON = 0;
if plotON
    



figure;
subplot(121);hold on;
theta_hat_true = @(x) x.^(-1/2);
plot(xgrid, real(theta_hat(xgrid)), 'DisplayName', 'est Lap[gamma]');
plot(xgrid, theta_hat_true(xgrid), 'DisplayName', 'True Lap[gamma]');
legend;

subplot(122);hold on;
theta_hat_true = @(x) x.^(-1/2);
xgrid_log = 10.^(linspace(-5, 5, 500));
plot(log10(xgrid_log), real(theta_hat(xgrid_log)), 'DisplayName', 'est Lap[gamma]');
plot(log10(xgrid_log), theta_hat_true(xgrid_log), 'DisplayName', 'True Lap[gamma]');
legend;


figure; 
subplot(121);hold on;

gamma = @(x) x.^(-1/2)/(sqrt(pi));    
xgrid = 0:0.01:20;

plot(xgrid, real(f(xgrid)), '-.', 'MarkerSize', 10, 'displayname', 'est')
plot(xgrid, gamma(xgrid), 'displayname', 'true Lap[gamma]')
legend;
ylim([0, 4])
title('Memory kernel')

subplot(122);hold on;grid on;
plot(log10(xgrid_log), log10(real(f(xgrid_log))), '-.', 'MarkerSize', 10, 'displayname', 'est')
plot(log10(xgrid_log), log10(gamma(xgrid_log)), 'displayname', 'true Lap[gamma]')
legend;
ylim([-2, 4])
title('Memory kernel, log scale')
    
end
end
   