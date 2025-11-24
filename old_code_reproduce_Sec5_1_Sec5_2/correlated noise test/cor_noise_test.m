clc
close all
clear all
% rng(1)
% This is a test about generating correlated noise with given correlation
% function gamma

%%
a = 1;
R = @(x) exp(-abs(x)*a);
S = @(w) 1/(2*pi)*2*a./(w.^2 + a^2);

wu = 636.6192;
N = 3321;
dw = wu/N;
wgrid = (0:dw:wu-dw)';

M = 16384;
assert(M > 2*N, 'M has to be larger than 2N')
dt = 0.002;
tgrid = 0:dt:M*dt;


%% plot
% wgrid_plot = linspace(-20, 20, 1000);
% L = 20;
% dt = 0.01;
% tgrid_plot = -L:dt:L;
% figure;
% subplot(121)
% plot(tgrid_plot, R(tgrid_plot));
% xlim([-5, 5])
% subplot(122)
% plot(wgrid_plot, S(wgrid_plot));
% % xlim([0, 20])



%%
phi_seq = rand(N, 1)*2*pi;

A_seq = zeros(N, 1);
for n = 2:N
    A_seq(n) = sqrt(2*S((n-1)*dw)*dw);
end




f = zeros(M+1, 1);
for m = 1:M
%     phi_seq = rand(N, 1)*2*pi;
    cos_seq = cos(wgrid*(m-1)*dt + phi_seq);
    f(m) = sqrt(2)*sum(A_seq.*cos_seq);
end
% 
% 
%%
A_seq_hat = zeros(N, 1);
for n = 2:N
    A_seq_hat(n) = sqrt(2*S(wgrid(n))*dw);
end
A_seq_hat(2) = sqrt((2*S(dw)+ 4/3*S(0))*dw);
A_seq_hat(3) = sqrt((2*S(2*dw)- 1/3*S(0))*dw);



f1 = zeros(M+1, 1);

% phi_seq = rand(N, 1)*2*pi;


for m = 1:M+1
%     phi_seq = rand(N, 1)*2*pi;
%     cos_seq = ;
    f1(m) = sqrt(2)*sum(A_seq_hat.*cos(wgrid*(m-1)*dt + phi_seq));
end


figure;hold on;
plot(tgrid, f1);
plot(tgrid, f);
%%
[acf1,lags1] = autocorr(f1, M);
[acf,lags] = autocorr(f, M);

figure;hold on;
% plot(lags, acf)
plot(tgrid, R(tgrid))
plot(tgrid, acf)
plot(tgrid, acf1)

legend('True', 'std', 'Better')
%


%%


% A_seq = zeros(M, 1);
% for n = 1:N
%     A_seq(n) = sqrt(2*S(wgrid(n))*dw);
% end
% 
% phi_seq = zeros(M, 1);
% phi_seq(1:N) = rand(N, 1)*2*pi;
% 
% B_seq = sqrt(2)*A_seq.*exp(1i*phi_seq);
% e_seq = exp(1i*(0:(M-1))*dw)';
% 
% dt = 0.002;
% tgrid = 0:dt:(M-1)*dt;
% f = zeros(M, 1);
% 
% for m = 1:M
%     phi_seq = zeros(M, 1);
%     phi_seq(1:N) = rand(N, 1)*2*pi;
%     
%     B_seq = sqrt(2)*A_seq.*exp(1i*phi_seq);
%     e_seq = exp(1i*(0:(M-1))*dw)';
%     f(m) = real(sum(B_seq.*(e_seq).^((m-1)*dt)));
% end
% 
% figure;
% plot(tgrid, f)



%
% %%
% L = 20;
% dt = 0.01;
% tgrid = -L:dt:L;
% N = length(tgrid);
%
%
% a = 1;
% f = @(x) exp(-abs(x)*a);
% % f = @(x) cos(10*x);
%
% X = f(tgrid);
% Y = fftshift(fft((X)));
%
%
% bd = 1/dt*2*pi;
% fgrid = (1:N)/N*bd - 0.5*bd;
%
%
% figure;
% subplot(121);
% plot(tgrid, X)
% subplot(122);
% plot(fgrid, abs(Y))
% xlim([-5, 5])


%%
% theta = randn


% %%
% dx = .1;
%  x = 0:dx:50;
% fs = 1/dx;
%    y = f(x);
%    Y = dx*fft(y); %Continous FT and DFT differ by scale dx
% X = -fs/2:fs/length(x):fs/2 - fs/length(x);
%
% Y = fftshift(Y);
% plot(X,real(Y));