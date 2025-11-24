fig3 = figure;
t = tiledlayout(1, 4, "TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';


xrange = 25;

nexttile;hold on;grid on;
plot(tgrid, corr_func.h, 'LineWidth', 3)
plot(tgrid, Prony_h(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.h(1:VACF_obs_len), '-o', 'LineWidth', 2)
plot(-10, -10, '.', 'MarkerSize', 0.001, 'color','white')
xlim([0, xrange])
% ylim([-1.5, 1.5])
title("h");
h_error_lgd = ['$\|\widetilde{h} - h_{true}\|_{L^2(\rho)} = $', num2str(h_error)];
legend("True", "Estimated", "Observed", h_error_lgd, 'location', 'southeast', 'interpreter' ,'Latex')
xlabel('time t')



nexttile;hold on;grid on;
plot(tgrid, corr_func.f, 'LineWidth', 3);
plot(tgrid, result_f.h(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.f(1:VACF_obs_len), '-o', 'LineWidth', 2)
plot(-10, -10, '.', 'MarkerSize', 0.001, 'color','white')
xlim([0, xrange])
% ylim([-0.1, 0.1])
title("f");
f_error_lgd = ['$\|\widetilde{f} - f_{true}\|_{L^2(\rho)} = $', num2str(f_error)];
legend("True", "Estimated", "Observed", f_error_lgd, 'location', 'southeast', 'interpreter' ,'Latex')
xlabel('time t')


nexttile;hold on;grid on;
plot(tgrid, corr_func.g, 'LineWidth', 3);
plot(tgrid, result_g.g(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.g(1:VACF_obs_len), '-o', 'LineWidth', 2)
plot(-10, -10, '.', 'MarkerSize', 0.001, 'color','white')
xlim([0, xrange])
title("g");
g_error_lgd = ['$\|\widetilde{g} - g_{true}\|_{L^2(\rho)} = $', num2str(g_error)];
legend("True", "Estimated", "Observed", g_error_lgd, 'location', 'southeast', 'interpreter' ,'Latex')
xlabel('time t')


nexttile;hold on;grid on;
plot(tgrid, corr_func.dg, 'LineWidth', 3);
plot(tgrid, result_g.dg(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.dg(1:VACF_obs_len), '-o', 'LineWidth', 2)
plot(-10, -10, '.', 'MarkerSize', 0.001, 'color','white')
xlim([0, xrange])
title("g'");
dg_error_lgd = strcat("$\|\widetilde{g}' - g_{true}'\|_{L^2(\rho)} = $", num2str(dg_error));
legend(["True", "Estimated", "Observed", dg_error_lgd], 'location', 'southeast', 'interpreter' ,'Latex')
xlabel('time t')






set(gcf, 'position', [300, 300, 1300, 400])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',14)



% print(fig3, [fig_dir, 'h_g_est.pdf'], '-dpdf', '-r0')
%%
