function f = generate_correlated_noise(R, S, wu, N, M, varargin)
% This function generate noise with given auto correlation
% Corresponding to the Generalized Langevin Equation with memory kernel
% Given autocorr R and its Fourier transform S ()
% This code generate a sample trajectory f, with mean 0 and auto coor R.
% The majority of the structure follows from M. Shinozuka & G. Deodatis, 1991
% With a bette convergence small change provided by B. Hu & W. Schiehlen, 1997

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



%% Parameters

p = inputParser;

addRequired(p, 'R');
addRequired(p, 'S');
addRequired(p, 'wu');
addRequired(p, 'N');
addRequired(p, 'M');

addOptional(p, 'dt', pi/wu);
addOptional(p, 'plotON', 0);


parse(p, R, S, wu, N, M, varargin{:});

dt      = p.Results.dt;
plotON  = p.Results.plotON;


dw          = wu/N;
T0          = 2*pi/dw;
wgrid       = (0:dw:wu-dw)';
T           = M*dt;

assert(dw<1, 'dw needs to be smaller than 1');
assert(dt <=pi/wu, 'dt nees to be smaller than pi/wu');

tgrid = (0:dt:T-dt)';

%% Some sample settings
% a = 1;
% R = @(x) exp(-abs(x)*a);
% S = @(w) 1/(2*pi)*2*a./(w.^2 + a^2);
% wu = 100*pi;
% N = 2000;
% M = 2^15;
%% Another sample settings
% R = @(x) (1-3*x.^2)./((1 + x.^2).^3);
% S = @(w) 1/4*w.^2.*exp(-abs(w));
%
% wu = 4*pi;
% N = 128;
% M = 2560;
% dt = 0.25;



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


%%

if plotON
    plot_noise(R, S, tgrid, f, wu)

    % %%
    % figure;grid on;
    % subplot(2, 2, 1);
    % plot(wgrid, S(wgrid));
    % xlim([0, 20])
    % title('Spectral function S')
    % 
    % % direct approximation
    % RR = zeros(M, 1);
    % for t = 1:M
    %     RR(t) = sum(f(1:end-t+1).*f(t:end))*dt/T;
    % end
    % 
    % % % use autocorr function
    % % [acf, ~] = autocorr(f, M-1);
    % % acf = acf * var(f);
    % 
    % % use autocox function
    % temp = xcov(f, f, 'unbiased');
    % temp = temp(M:end);
    % 
    % 
    % subplot(2, 2, 2);hold on;
    % plot(tgrid, temp, 'DisplayName','XCOV');
    % plot(tgrid, R(tgrid), 'LineWidth', 4, 'DisplayName','True');
    % % plot(tgrid, acf, '.', 'markersize', 15, 'DisplayName','ACF');
    % plot(tgrid, RR, '.', 'markersize', 15, 'DisplayName','direct');
    % 
    % 
    % xlim([0, 20])
    % title('Autocorrelation')
    % legend()
    % % legend('true', 'approx');
    % 
    % 
    % 
    % subplot(2, 2, [3, 4]);
    % hold on;
    % plot(tgrid, f);
    % xlabel('t')
    % title('Trajectory');

end

end






