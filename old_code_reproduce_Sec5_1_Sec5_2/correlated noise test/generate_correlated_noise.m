function f = generate_correlated_noise(R, S, wu, N, M, dt)
% This function generate noise with given auto correlation
% Corresponding to the Generalized Langevin Equation with memory kernel
% Given autocorr R and its Fourier transform S ()
% This code generate a sample trajectory f, with mean 0 and auto coor R.
% The majority of the structure follows from M. Shinozuka & G. Deodatis, 1991
% With a better convergence small change provided by B. Hu & W. Schiehlen, 1997

% Input:    R: auto correlation function
%           S: Fourier transform of R
%           Wu: Upper bound for the spectral integral
%           N: Number of grid points in the spectral domain
%           M: Number of grid points in the time domain

% Output: f: Stochastic process with mean 0 and autocorr R.

% Requirement:  1. The result f has period T0 = 2*pi/dw, dw = Wu/N
%               2. dt needs to be greater than pi/Wu
%               3. Only when dw < 1, the better convergence will work
%               4. When M*dt = T == T0, or T --> inf,
%                   the autocorr of the generated process has a better fit.

% Plot options
plotON = 1;
% Plot R, S, f, and autocorr computed from f.

%% Some sample settings
% a = 1;
% R = @(x) exp(-abs(x)*a);
% S = @(w) 1/(2*pi)*2*a./(w.^2 + a^2);
% wu = 100*pi;
% N = 2000;
% M = 2^15;
% dt = 0.01;
%% Another sample settings
% R = @(x) (1-3*x.^2)./((1 + x.^2).^3);
% S = @(w) 1/4*w.^2.*exp(-abs(w));
% 
% wu = 4*pi;
% N = 128;
% M = 2560;
% dt = 0.25;

%% Parameters
if nargin == 5
    dt = pi/wu;     % Usually we set dt = pi/Wu. dt can be specified in the command.
end

dw = wu/N;
wgrid = (0:dw:wu-dw)';
T0 = 2*pi/dw;


T = M*dt;

assert(dw<1, 'dw needs to be smaller than 1');
assert(dt <=pi/wu, 'dt nees to be smaller than pi/wu');

tgrid = 0:dt:T-dt;


%%
% random noise
phi_seq = rand(N, 1)*2*pi;

% spectral function evaluated on wgrid
A_seq = zeros(N, 1);
for n = 2:N
    A_seq(n) = sqrt(2*S((n-1)*dw)*dw);
end

% Better convergence, introduced in B. Hu & W. Schiehlen, 1997
A_seq(2) = sqrt((2*S(dw)+ 4/3*S(eps))*dw);
A_seq(3) = sqrt((2*S(2*dw)- 1/3*S(eps))*dw);

% FFT speed up
B_seq = sqrt(2)*A_seq.*exp(1i*phi_seq);
f = zeros(M, 1);

% spectral integration
for p = 1:M
    E = exp(1i*wgrid*(p-1)*dt);
    I = sum(B_seq.*E);
    f(p) = real(I);
end



if plotON
    %%
    figure;
    subplot(2, 2, 1);
    plot(wgrid, S(wgrid));
    xlim([0, 20])
    title('Spectral function S')
    
    subplot(2, 2, 2);hold on;
    RR = zeros(M, 1);
    for t = 1:M
        RR(t) = sum(f(1:end-t+1).*f(t:end))*dt/T;
    end
%     [acf, ~] = autocorr(f, M-1);
    plot(tgrid, R(tgrid), 'LineWidth', 4);
%     plot(tgrid, acf, '.', 'markersize', 15);
    plot(tgrid, RR, '.', 'markersize', 15);
    
    
    xlim([0, 8])
    title('Autocorrelation')
    legend('true', 'approx');
    
    
    
    subplot(2, 2, [3, 4]);
    hold on;
    plot(tgrid, f);
    xlabel('t')
    title('Trajectory');
    
    
end

end