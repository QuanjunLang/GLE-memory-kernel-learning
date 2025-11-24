% p = 5;
% w_seq = 3*(rand(p, 1));
% lambda_seq = rand(p, 1);

p = 1;
q = 3;

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







R = @(x) 0*x;
S = @(x) 0*x;
for k = 1:2*p+q
    R = @(x) R(x) + w_seq(k)*exp(lambda_seq(k)*x);
    S = @(x) S(x) + w_seq(k)*lambda_seq(k)/pi ./(x.^2 + lambda_seq(k)^2);
end

a = 1;
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


figure;
subplot(121)
fplot(R, [0, 15], 'LineWidth', 4);title('True Kernel \gamma')
subplot(122)
fplot(R_lap, [0, 15], 'LineWidth', 4);title('Laplace transform of \gamma')

a = 1;