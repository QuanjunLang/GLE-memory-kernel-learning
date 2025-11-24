clc
close all
clear all




%%
a = 1;
R = @(x) (1-3*x.^2)./((1 + x.^2).^3);
S = @(w) 1/4*w.^2.*exp(-abs(w));

wu = 4*pi;
N = 128;
dw = wu/N;
wgrid = (0:dw:wu-dw)';

% M = 2^13;
% assert(M > 2*N, 'M has to be larger than 2N')
% dt = 0.001;
% tgrid = 0:dt:(M-1)*dt;
% T = M*dt;


T0 = 2*pi/dw;
% T = T0;
M = 2560;
dt = 0.25;
T = M*dt;
tgrid = 0:dt:T-dt;
% M = M+1;
% dt = tgrid(2) - tgrid(1);

pi/wu
dt

%%
phi_seq = rand(N, 1)*2*pi;

A_seq = zeros(N, 1);
for n = 2:N
    A_seq(n) = sqrt(2*S((n-1)*dw)*dw);
end
% A_seq(2) = sqrt((2*S(dw)+ 4/3*S(0))*dw);
% A_seq(3) = sqrt((2*S(2*dw)- 1/3*S(0))*dw);

B_seq = sqrt(2)*A_seq.*exp(1i*phi_seq);
% B_seq(N+1:end) = 0;


f_FFT = zeros(M, 1);
f_STD = zeros(M, 1);
for p = 1:M
    E = exp(1i*wgrid*(p-1)*dt);
    I = sum(B_seq.*E);
    f_FFT(p) = real(I);
    
    cos_seq = cos(wgrid*(p-1)*dt + phi_seq);
    f_STD(p) = sqrt(2)*sum(A_seq.*cos_seq);
end

figure;hold on;
plot(f_FFT);
plot(f_STD, '.', 'LineWidth', 20);

%%
RR = zeros(M, 1);
for t = 1:M 
    RR(t) = sum(f_STD(1:end-t+1).*f_STD(t:end))*dt/T;
end



%%
[acf,lags] = autocorr(f_FFT, M-1);

figure;hold on;
plot(tgrid, acf)
plot(tgrid, R(tgrid))
plot(tgrid, RR, '.', 'LineWidth' ,20)
xlim([0, 5])


