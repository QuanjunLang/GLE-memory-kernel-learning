fig3 = figure;
t = tiledlayout(2, 3, "TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';



nexttile;hold on;grid on;
plot(tgrid, corr_func.h, 'LineWidth', 3)
plot(tgrid, Prony_h(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.h(1:VACF_obs_len), '-o', 'LineWidth', 2)
xlim([0, 20])
title("h: Autocorrelation of v");
legend('True h', 'Estimated h', 'Observed h')
xlabel('time t')

nexttile;hold on;grid on;
plot(tgrid, corr_func.dh, 'LineWidth', 3);
plot(tgrid, result_h.dh(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.dh(1:VACF_obs_len), '-o', 'LineWidth', 2)
xlim([0, 20])
title("h'");
legend("True g'", "Estimated g'", "Observed g'", 'location', 'southeast')
xlabel('time t')


nexttile;hold on;grid on;
plot(tgrid, corr_func.ddh, 'LineWidth', 3);
plot(tgrid, result_h.ddh(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.ddh(1:VACF_obs_len), '-o', 'LineWidth', 2)
xlim([0, 20])
title("h''");
legend("True g'", "Estimated g'", "Observed g'", 'location', 'southeast')
xlabel('time t')


nexttile;hold on;grid on;
plot(tgrid, corr_func.g, 'LineWidth', 3);
plot(tgrid, result_g.g(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.g(1:VACF_obs_len), '-o', 'LineWidth', 2)
xlim([0, 20])
title("g");
legend("True", "Estimated ", "Observed ", 'location', 'southeast')
xlabel('time t')

nexttile;hold on;grid on;
plot(tgrid, corr_func.dg, 'LineWidth', 3);
plot(tgrid, result_g.dg(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.dg(1:VACF_obs_len), '-o', 'LineWidth', 2)
xlim([0, 20])
title("g'");
legend("True g'", "Estimated g'", "Observed g'", 'location', 'southeast')
xlabel('time t')


nexttile;hold on;grid on;
plot(tgrid, corr_func.f, 'LineWidth', 3);
plot(tgrid, result_f.h(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.f(1:VACF_obs_len), '-o', 'LineWidth', 2)
xlim([0, 20])
title("esternal forcing term");
legend("True f", "Estimated f", "Observed f", 'location', 'southeast')
xlabel('time t')



set(gcf, 'position', [300, 300, 800, 200])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',14)


print(fig3, [fig_dir, 'h_g_est.pdf'], '-dpdf', '-r0')