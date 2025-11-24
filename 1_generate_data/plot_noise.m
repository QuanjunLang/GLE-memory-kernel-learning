function plot_noise(R, S, tgrid, f, wu)
% plot the noise and the autocorrelation

M = length(f);
dt = tgrid(2) - tgrid(1);
T = tgrid(end);
%%
figure;grid on;
subplot(2, 2, 1);
fplot(S);
xlim([0, wu])
title('Spectral function S')

% direct approximation
RR = zeros(M, 1);
for t = 1:M
    RR(t) = sum(f(1:end-t+1).*f(t:end))*dt/T;
end




% use autocox function
temp = xcov(f, f, 'unbiased');
temp = temp(M:end);


subplot(2, 2, 2);hold on;grid on;
idx = find(abs(R(tgrid)) > 1e-3, 1, 'last');
plot(tgrid, R(tgrid), 'LineWidth', 4, 'DisplayName','True');
plot(tgrid, temp, 'DisplayName','XCOV');
plot(tgrid, RR, '.', 'markersize', 15, 'DisplayName','direct');


% % use autocorr function
% [acf, ~] = autocorr(f, M-1);
% acf = acf * var(f);
% plot(tgrid, acf, '.', 'markersize', 15, 'DisplayName','ACF');

xlim([0, tgrid(idx)])
title('Autocorrelation')
legend()


subplot(2, 2, [3, 4]);
hold on;
plot(tgrid, f);
xlabel('t')
title('Trajectory');

end