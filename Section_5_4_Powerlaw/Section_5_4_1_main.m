
close all
clear all

result_exp_L_216_M_1            = Test_exp_L_216_M_1();
result_exp_L_212_M_2000         = Test_exp_L_212_M_2000();
result_powerlaw_L_216_M_1       = Test_powerlaw_L_216_M_1();
result_powerlaw_L_212_M_2000    = Test_powerlaw_L_212_M_2000();


%%
tgrid = linspace(0, 10, 500)';  % time grid


figure;
t = tiledlayout(1, 2, "TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';

lgd216 = sprintf('$L = 2^{16}, M = 1,  \\quad  \\ \\, \\|\\widehat{h} - h\\|_{L^2(\\rho)} = %.4f$', result_exp_L_216_M_1.h_error);
lgd212 = sprintf('$L = 2^{12}, M = 2000, \\|\\widehat{h} - h\\|_{L^2(\\rho)} = %.4f$', result_exp_L_212_M_2000.h_error );


nexttile;hold on;grid on;
plot(tgrid, result_exp_L_216_M_1.h_true(tgrid), '-' , 'LineWidth', 2, 'DisplayName', 'True $h$');
plot(tgrid, result_exp_L_216_M_1.h_est(tgrid),      '-.', 'LineWidth', 4, 'DisplayName', lgd216);
plot(tgrid, result_exp_L_212_M_2000.h_est(tgrid),   '-.', 'LineWidth', 4, 'DisplayName', lgd212)
xlim([0, min(20, tgrid(end))])
title("h(t)")
legend('Interpreter','latex')


lgd216 = sprintf('$L = 2^{16}, M = 1,  \\quad  \\ \\, \\|\\widehat{g} - g\\|_{H^1_\\alpha(\\rho)} = %.4f$', result_exp_L_216_M_1.g_error);
lgd212 = sprintf('$L = 2^{12}, M = 2000, \\|\\widehat{g} - g\\|_{H^1_\\alpha(\\rho)} = %.4f$', result_exp_L_212_M_2000.g_error );

nexttile;hold on;grid on;
plot(tgrid, result_exp_L_216_M_1.g_true(tgrid), '-' , 'LineWidth', 2, 'DisplayName', 'True $g$')
plot(tgrid, result_exp_L_216_M_1.g_est(tgrid),      '-.', 'LineWidth', 4, 'DisplayName', lgd216);
plot(tgrid, result_exp_L_212_M_2000.g_est(tgrid),   '-.', 'LineWidth', 4, 'DisplayName', lgd212);
xlim([0, min(20, tgrid(end))])
title("g(t)")
legend('Interpreter','latex')

% 
% nexttile;hold on;grid on;
% plot(tgrid, result_exp_L_216_M_1.dg_true(tgrid), '-' , 'LineWidth', 3, 'DisplayName', "True $g'$")
% plot(tgrid, result_exp_L_216_M_1.dg_est(tgrid),      '-.', 'LineWidth', 4, 'DisplayName', '$L = 2^{16}, M = 1$')
% plot(tgrid, result_exp_L_212_M_2000.dg_est(tgrid),   '-.', 'LineWidth', 4, 'DisplayName', '$L = 2^{12}, M = 2000$')
% xlim([0, min(20, tgrid(end))])
% title("$g'(t)$", 'Interpreter','latex')
% legend('Interpreter','latex')


set(gcf, 'position', [300, 300, 1000, 220])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',14)

saveas(gcf, 'exp_vacf.pdf')

%%
figure;
t = tiledlayout(1, 2, "TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';

rho = result_exp_L_216_M_1.rho;

nexttile;hold on;
plot_temp(tgrid, rho, result_exp_L_216_M_1, 'Estimation result, L = 2^{16}, M = 1')

nexttile;hold on;
plot_temp(tgrid, rho, result_exp_L_212_M_2000, 'Estimation result, L = 2^{12}, M = 2000')


set(gcf, 'position', [300, 300, 1000, 300])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',15)

saveas(gcf, 'exp_result.pdf')
%%
tgrid = linspace(0, 10, 500)';  % time grid

figure;
t = tiledlayout(1, 2, "TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';


nexttile;hold on;grid on;
% plot(tgrid, result_powerlaw_L_216_M_1.h_true(tgrid), '-' , 'LineWidth', 3, 'DisplayName', 'True $h$');
plot(tgrid, result_powerlaw_L_216_M_1.h_est(tgrid),      '-.', 'LineWidth', 4, 'color', [0.8500, 0.3250, 0.0980], 'DisplayName', '$L = 2^{16}, M = 1$')
plot(tgrid, result_powerlaw_L_212_M_2000.h_est(tgrid),   '-', 'LineWidth', 2, 'color', [0.9290 0.6940 0.1250], 'DisplayName', '$L = 2^{12}, M = 2000$')
xlim([0, min(20, tgrid(end))])
title('h(t)')
legend('Interpreter','latex')


nexttile;hold on;grid on;
% plot(tgrid, result_powerlaw_L_216_M_1.g_true(tgrid), '-' , 'LineWidth', 3, 'DisplayName', 'True $g$')
plot(tgrid, result_powerlaw_L_216_M_1.g_est(tgrid),      '-.', 'LineWidth', 4, 'color', [0.8500, 0.3250, 0.0980], 'DisplayName', '$L = 2^{16}, M = 1$')
plot(tgrid, result_powerlaw_L_212_M_2000.g_est(tgrid),   '-', 'LineWidth', 2,'color', [0.9290 0.6940 0.1250], 'DisplayName','$L = 2^{12}, M = 2000$')
xlim([0, min(20, tgrid(end))])
title("g(t)")
legend('Interpreter','latex')


% nexttile;hold on;grid on;
% % plot(tgrid, result_powerlaw_L_216_M_1.dg_true(tgrid), '-' , 'LineWidth', 3, 'DisplayName', "True $g'$")
% plot(tgrid, result_powerlaw_L_216_M_1.dg_est(tgrid),      '-.', 'LineWidth', 4, 'DisplayName', '$L = 2^{16}, M = 1$')
% plot(tgrid, result_powerlaw_L_212_M_2000.dg_est(tgrid),   '-.', 'LineWidth', 4, 'DisplayName', '$L = 2^{12}, M = 2000$')
% xlim([0, min(20, tgrid(end))])
% title("$g'(t)$", 'Interpreter','latex')
% legend('Interpreter','latex')


set(gcf, 'position', [300, 300, 1000, 220])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',14)

saveas(gcf, 'powerlaw_vacf.pdf')
%%
figure;
t = tiledlayout(1, 2, "TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';

rho = result_exp_L_216_M_1.rho;

nexttile;hold on;
plot_temp(tgrid, rho, result_powerlaw_L_216_M_1, 'Estimation result, L = 2^{16}, M = 1')

nexttile;hold on;
plot_temp(tgrid, rho, result_powerlaw_L_212_M_2000, 'Estimation result, L = 2^{12}, M = 2000')


set(gcf, 'position', [300, 300, 1000, 300])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',15)




saveas(gcf, 'powerlaw_result.pdf')