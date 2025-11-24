clc
close all
clear all
addPathAll;

result_L_216_M_1        = test_L_216_M_1();
result_L_212_M_1        = test_L_212_M_1();
result_L_212_M_2000     = test_L_212_M_2000();


%%
close all

omega = result_L_216_M_1.omega;
rho = @(x) exp(-2*omega*x);

theta_L212_M_1      = result_L_212_M_1.theta;
theta_L216_M_1      = result_L_216_M_1.theta;
theta_L212_M_2000   = result_L_212_M_2000.theta;

theta_true          = result_L_212_M_1.theta_true;

lgd_L212_M_1 = sprintf('$L = 2^{12}, M = 1, \\quad  \\ \\, \\|\\theta - \\gamma\\|_{L^2(\\rho)} = %.4f$', result_L_212_M_1.theta_error);
lgd_L216_M_1 = sprintf('$L = 2^{16}, M = 1, \\quad  \\ \\, \\|\\theta - \\gamma\\|_{L^2(\\rho)} = %.4f$', result_L_216_M_1.theta_error);
lgd_L212_M_2000 = sprintf('$L = 2^{12}, M = 2000, \\|\\theta - \\gamma\\|_{L^2(\\rho)} = %.4f$', result_L_212_M_2000.theta_error);

figure;
hold on;
left_color = [.5 .5 0];
right_color = [0 .5 .5];
set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
grey = 0.7*ones(1, 3);

tgrid = linspace(0, 16, 500)';  % time grid





yyaxis right
fill_rho_xgrid = [tgrid; fliplr(tgrid')'];
fill_rho_val = [rho(tgrid); zeros(size(tgrid))];
fill(fill_rho_xgrid, fill_rho_val, 'k', 'FaceAlpha', 0.15, 'DisplayName','Measure $\rho$', 'EdgeColor', 'none')
ylim([0, 2])
ylabel('measure \rho')



yyaxis left
grid on;
plot(tgrid, theta_true(tgrid),          '-', 'LineWidth', 5, 'Color', [0 0.45 0.74], 'DisplayName', 'True kernel $\gamma$')
plot(tgrid, theta_L212_M_1(tgrid),      '-.', 'LineWidth', 5, 'color', [0.4940    0.1840    0.5560] , 'DisplayName', lgd_L212_M_1);
plot(tgrid, theta_L216_M_1(tgrid),      '-',  'LineWidth', 5, 'color', [0.4660    0.6740    0.1880], 'DisplayName', lgd_L216_M_1);
plot(tgrid, theta_L212_M_2000(tgrid),   '-.',  'LineWidth', 5, 'color', [0.8500    0.3250 0.0980], 'DisplayName', lgd_L212_M_2000);

ylabel('memory kernel')



ylim([min(theta_true(tgrid))*1.5, max(theta_true(tgrid))*1.5])
legend('Interpreter','latex')

title('Estimation result of the memory kernel')


ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = grey;

set(gcf, 'position', [300, 300, 700, 300])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',15)



%%
saveas(gcf, 'External_force_531.pdf');  % gcf = get current figure handle