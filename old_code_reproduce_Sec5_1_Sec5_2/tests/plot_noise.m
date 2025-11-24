function fig = plot_noise(R, S, M, para, f, xrange)

wu = para.wu;
N = para.N;
dt = pi/wu;
dw = wu/N;
wgrid = (0:dw:wu-dw)';
T = M*dt;
tgrid = 0:dt:T-dt;
% f = f * para.beta;
R = @(x) R(x)/(para.beta^2);
%%
fig = figure;
subplot(2, 2, 2);
plot(wgrid, S(wgrid), 'LineWidth', 4);
xlim([0, 10])
title('Spectral function S')

subplot(2, 2, 1);hold on;
RR = zeros(M, 1);
for t = 1:M/10
    RR(t) = sum(f(1:end-t+1).*f(t:end))*dt/T;
end
%     [acf, ~] = autocorr(f, M-1);
plot(tgrid, R(tgrid), 'LineWidth', 4);
%     plot(tgrid, acf, '.', 'markersize', 15);
% plot(tgrid(1:3:end), RR(1:3:end), 'o', 'markersize', 10);
plot(tgrid(1:5:end), RR(1:5:end), '-.', 'LineWidth', 4, 'markersize', 3);

xlim([0, xrange])
title('Autocorrelation')
legend('true', 'approximate');



subplot(2, 2, [3, 4]);
hold on;
plot(tgrid, f);
xlabel('t')
xlim([0, T])
title('Correlated Noise');


set(findall(gcf,'-property','FontSize'),'FontSize',12)
