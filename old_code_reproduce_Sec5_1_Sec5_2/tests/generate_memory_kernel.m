function [R, S, R_lap, w_seq, lambda_seq, h_true, h_true_lap, dh_true, ddh_true] = generate_memory_kernel(para)

p = para.p;
q = para.q;
% plotON = para.plotON;

% p = 1;
% q = 3;

w_seq = zeros(2*p+q, 1);
lambda_seq = zeros(2*p+q, 1);

for k = 1:p
    w_seq(2*k-1) = (rand()-0.5);
    lambda_seq(2*k-1) = (rand()-0.5)*1i - 3*rand();
    w_seq(2*k) = w_seq(2*k-1);
    lambda_seq(2*k) = conj(lambda_seq(2*k-1));
end

for k = 2*p+1:2*p+q
    w_seq(k) = rand();
    lambda_seq(k) = -3*rand();
end

% beta = 5;
% w_seq = w_seq/beta;




R = @(x) 0*x;
S = @(x) 0*x;
for k = 1:2*p+q
    R = @(x) R(x) + w_seq(k)*exp(lambda_seq(k)*x);
    S = @(x) S(x) + w_seq(k)*(-lambda_seq(k))/pi ./(x.^2 + lambda_seq(k)^2);
end

% a = 1;
sgn = sign(R(0));
R = @(x) R(x)*sgn;
S = @(x) S(x)*sgn;
% S is the Fourier transform of R
% The fourier transform used here is 
% F[R] = 1/2pi int R(x) exp(-iwx)dx = S(w)

R_lap = @(x) 0*x;

for k = 1:2*p+q
    R_lap = @(x) R_lap(x) + w_seq(k)./(x - lambda_seq(k));
end

R_lap = @(x) R_lap(x)*sgn;


if para.plotON
range = 30;
figure;
subplot(121)
fplot(R, [0, range], 'LineWidth', 4);title('True Kernel \gamma')
subplot(122)
fplot(R_lap, [0, range], 'LineWidth', 4);title('Laplace transform of \gamma')

end

%% getting true h, when F = 0
syms s
k = length(w_seq);

R_lap_sym = 0*s;
for i = 1:k
    R_lap_sym = R_lap_sym + w_seq(i)/(s - lambda_seq(i));
end

h_lap_sym = vpa(1/(s + R_lap_sym));      
h_lap_sym = partfrac(h_lap_sym, 'FactorMode', 'complex');
% pretty(h_lap_sym)
C = children(h_lap_sym);

for i = 1:length(C)
    [N, D] = numden(C{i});
    coef = coeffs(D);
    par = coef(2);
    h_w_seq(i) = N/par;
    h_lam_seq(i) = coef(1)/par;
end
h_lam_seq = -h_lam_seq;
% r = exp(h_lam_seq);


h_true = 0*s;
for i = 1:length(C)
    h_true = h_true + h_w_seq(i)*exp(s*h_lam_seq(i));
end
h_true = matlabFunction(h_true);
h_true = @(x) real(h_true(x));

dh_true = 0*s;
for i = 1:length(C)
    dh_true = dh_true + h_lam_seq(i)*h_w_seq(i)*exp(s*h_lam_seq(i));
end
dh_true = matlabFunction(dh_true);
dh_true = @(x) real(dh_true(x));

ddh_true = 0*s;
for i = 1:length(C)
    ddh_true = ddh_true + h_lam_seq(i)*h_lam_seq(i)*h_w_seq(i)*exp(s*h_lam_seq(i));
end
ddh_true = matlabFunction(ddh_true);
ddh_true = @(x) real(ddh_true(x));

%%
h_true_lap = 0*s;
for i = 1:length(C)
    h_true_lap = h_true_lap + h_w_seq(i)/(s - h_lam_seq(i));
end

h_true_lap = matlabFunction(h_true_lap);
h_true_lap = @(x) real(h_true_lap(x));

end