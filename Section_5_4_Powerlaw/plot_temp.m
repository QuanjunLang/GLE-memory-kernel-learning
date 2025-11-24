function plot_temp(tgrid, rho, Info, ttl)

colors.blue     = [0, 0.4470, 0.7410];
colors.red      = [0.8500, 0.3250, 0.0980];
colors.yellow   = [0.9290, 0.6940, 0.1250];
colors.purple   = [0.4940, 0.1840, 0.5560];
colors.green    = [0.4660, 0.6740, 0.1880];
colors.cyan     = [0.3010, 0.7450, 0.9330];
colors.maroon   = [0.6350, 0.0780, 0.1840];
colors.purple   = [0.4940    0.1840    0.5560];
colors.yellow   = [0.9290 0.6940 0.1250];


hold on;

theta_0      = Info.theta_0;
theta_1      = Info.theta_1;
theta_2      = Info.theta_2;
theta_L      = Info.theta_L;
theta_F      = Info.theta_F;
theta_T      = Info.theta_T;


lgd_T = 'True $\gamma$';
lgd_0 = sprintf('$\\theta_0, \\quad  \\ \\, \\|\\theta_0 - \\gamma\\|_{L^2(\\rho)} = %.4f$', Info.error_0 );
lgd_1 = sprintf('$\\theta_1, \\quad  \\ \\, \\|\\theta_1 - \\gamma\\|_{L^2(\\rho)} = %.4f$', Info.error_1 );
lgd_2 = sprintf('$\\theta_2, \\quad  \\ \\, \\|\\theta_2 - \\gamma\\|_{L^2(\\rho)} = %.4f$', Info.error_2 );
lgd_L = sprintf('$\\theta_L, \\quad  \\ \\, \\|\\theta_L - \\gamma\\|_{L^2(\\rho)} = %.4f$', Info.error_L );
lgd_F = sprintf('$\\theta_F, \\quad  \\ \\, \\|\\theta_F - \\gamma\\|_{L^2(\\rho)} = %.4f$', Info.error_F );


left_color = [.5 .5 0];
right_color = [0 .5 .5];
set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
grey = 0.7*ones(1, 3);


yyaxis right
fill_rho_xgrid = [tgrid; fliplr(tgrid')'];
fill_rho_val = [rho(tgrid); zeros(size(tgrid))];
fill(fill_rho_xgrid, fill_rho_val, 'k', 'FaceAlpha', 0.15, 'DisplayName','Measure $\rho$', 'EdgeColor', 'none')

% ylabel('measure \rho')



yyaxis left
grid on;
plot(tgrid, theta_T(tgrid),      '-', 'LineWidth', 5, 'Color',  colors.blue                 , 'DisplayName', lgd_T)
plot(tgrid, theta_0(tgrid),      '-.', 'LineWidth', 5, 'color', colors.red , 'DisplayName', lgd_0);
plot(tgrid, theta_1(tgrid),      '-', 'LineWidth', 5, 'color', colors.green , 'DisplayName', lgd_1);
plot(tgrid, theta_2(tgrid),      '-',  'LineWidth', 3, 'color', colors.purple, 'DisplayName', lgd_2);
plot(tgrid, theta_L(tgrid),      '.', 'MarkerSize', 10, 'linewidth', 4, 'color', [0.9290 0.6940 0.1250], 'DisplayName',lgd_L);
plot(tgrid, theta_F(tgrid),      '-.',  'LineWidth', 3, 'color', colors.cyan, 'DisplayName', lgd_F);

% ylabel('memory kernel')



ylim([min(theta_T(tgrid))-0.25, max(theta_T(tgrid))*1.5])
legend('Interpreter','latex')

title(ttl)


ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = grey;


end