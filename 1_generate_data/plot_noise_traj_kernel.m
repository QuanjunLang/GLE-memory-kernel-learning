function [fig, fig2] = plot_noise_traj_kernel(R, S, M, para, f, xrange, v, v_obs, tgrid_obs, T0, fig_dir)

wu = para.wu;
N = para.N;
dt = pi/wu;
dw = wu/N;
wgrid = (0:dw:wu-dw)';
T = M*dt;
tgrid = 0:dt:T-dt;
% f = f * para.beta;
% R = @(x) R(x)/(para.beta^2);

RR = zeros(M, 1);
for t = 1:M/10
    RR(t) = sum(f(1:end-t+1).*f(t:end))*dt/T;
end


% %%
% fig = figure;
% subplot(2, 2, 2);
% plot(wgrid, S(wgrid), 'LineWidth', 4);
% xlim([0, 10])
% title('Spectral function S')
% 
% subplot(2, 2, 1);hold on;
% 
% %     [acf, ~] = autocorr(f, M-1);
% plot(tgrid, R(tgrid), 'LineWidth', 4);
% %     plot(tgrid, acf, '.', 'markersize', 15);
% % plot(tgrid(1:3:end), RR(1:3:end), 'o', 'markersize', 10);
% plot(tgrid(1:5:end), RR(1:5:end), '-.', 'LineWidth', 4, 'markersize', 3);
% 
% xlim([0, xrange])
% title('Autocorrelation')
% legend('true', 'approximate');
% 
% 
% 
% subplot(2, 2, [3, 4]);
% hold on;
% plot(tgrid, f);
% xlabel('t')
% xlim([0, T])
% title('Correlated Noise');
% 
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% 
% 
% figure;hold on;title('Trajectory v')
% plot(tgrid, v, 'linewidth', 1);
% plot(tgrid_obs, v_obs, '.', 'markersize', 4);
% xline(T0, 'linewidth', 3)
% xlabel('time t')
% legend('True', 'obs', 'Period T0')
% 
% a = 1;
%% 
% fig = figure;
% t = tiledlayout(2, 3,"TileSpacing","compact");
% t.TileSpacing = 'compact';
% t.Padding = 'none';
% nexttile;hold on;
% grid on ;
% plot(tgrid, R(tgrid), 'LineWidth', 3);
% plot(tgrid(1:5:end), RR(1:5:end), '-.', 'LineWidth', 4, 'markersize', 3);
% xlim([0, xrange])
% title('Memory kernel')
% legend('True kernel \gamma', 'Autocorr of noise', 'location', 'northeast');
% 
% nexttile(2, [1, 2])
% 
% hold on;title('Trajectory v')
% plot(tgrid, v, 'linewidth', 1);
% plot(tgrid_obs, v_obs, '.', 'markersize', 4);
% xlim([0, T])
% 
% 
% legend('True', 'Observed')
% 
% nexttile;hold on;
% grid on;
% plot(wgrid, S(wgrid), 'LineWidth', 3);
% xlim([0, 10])
% legend('$\mathcal{F}(\gamma)$', 'location', 'southeast', 'interpreter', 'latex')
% title('Spectral function')
% 
% nexttile(5, [1,2])
% hold on;
% plot(tgrid, f, 'linewidth', 0.8);
% % xlabel('t')
% xlim([0, T])
% title('Correlated Noise');
% xline(T0, 'linewidth', 3)
% legend('Noise', 'Period T_0');
% xlabel('time t')
% set(gcf, 'PaperSize',[9, 4.5])
% set(gcf, 'position', [300, 300, 550, 300])
% set(findall(gcf,'-property','FontSize'),'FontSize',12)


%%
fig = figure;
t = tiledlayout(1, 2,"TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';




nexttile;hold on;
grid on ;
plot(tgrid, R(tgrid), 'LineWidth', 3);
plot(tgrid(1:5:end), RR(1:5:end), '-.', 'LineWidth', 4, 'markersize', 3);
xlim([0, xrange])
title('Memory kernel')
legend('True kernel \gamma', 'Auto cov of noise', 'location', 'northeast');
xlabel('time t')


nexttile;hold on;
grid on;
plot(wgrid, S(wgrid), 'LineWidth', 3);
xlim([0, 10])
xlabel('')
% legend('$\mathcal{F}[\gamma]$', 'location', 'northeast', 'interpreter', 'latex')
title('Spectral function')
xlabel('frequency')

set(fig, 'position', [300, 300, 400, 200])
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(fig,'-property','FontSize'),'FontSize', 14)








fig2 = figure;
t = tiledlayout(1, 4, "TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';
nexttile(3, [1, 2])

hold on;title('Trajectory v')
plot(tgrid, v, 'linewidth', 1);
plot(tgrid_obs, v_obs, '.', 'markersize', 8, 'color', [0.9290 0.6940 0.1250]);
xlim([0, T])
xlabel('time t')

legend('True', 'Observed')



nexttile(1, [1, 2])
hold on;
plot(tgrid, f, 'linewidth', 0.8);
% xlabel('t')
xlim([0, T])
title('Correlated Noise R');
xline(T0, 'linewidth', 3)
legend('Noise', 'Period T_0');
xlabel('time t')
% set(gcf, 'PaperSize',[9, 4.5], 'PaperUnits', 'normalized')
set(fig2, 'position', [300, 300, 800, 200])
set(fig2,'Units','Inches');
pos = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

set(findall(gcf,'-property','FontSize'),'FontSize',14)



end