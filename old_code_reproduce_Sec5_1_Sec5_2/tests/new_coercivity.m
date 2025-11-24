%% Coercivity analysis

% omega = 0.0001
% omega = 0.0000001;
% omega = 0.1;
tau_grid_real = 10.^linspace(-10, 10, 10000);
tau_grid = omega + 1i*tau_grid_real;

z_gamma_z = abs(tau_grid.*R_lap(tau_grid));
gamma_z = abs(R_lap(tau_grid));

gamma_z_comb = comb_coef*gamma_z.^2 + (1-comb_coef)*z_gamma_z.^2;
m_gamma = min(gamma_z_comb);
M_gamma = max(gamma_z_comb);

% Secondly, the bound for laplace transform of h_eps
% tau_grid = omega + 1i*tau_grid_real;
z_h_eps_z = abs(tau_grid.*result_h.h_lap(tau_grid));
h_eps_z = abs(result_h.h_lap(tau_grid));

m_z_h_eps_z = min(z_h_eps_z.^2);
M_z_h_eps_z = max(z_h_eps_z.^2);
m_h_eps_z = min(h_eps_z.^2);
M_h_eps_z = max(h_eps_z.^2);

h_eps_z_comb = comb_coef*h_eps_z.^2 + (1-comb_coef)*z_h_eps_z.^2;
m_h_eps_z_comb = min(h_eps_z_comb);
M_h_eps_z_comb = max(h_eps_z_comb);


fig_coer = figure;
t = tiledlayout(1, 4, "TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';


nexttile
hold on;grid on;
plot(log10(tau_grid_real), log10(gamma_z_comb), 'LineWidth', 4)
yline(log10(M_gamma), 'LineWidth', 4)
yline(log10(m_gamma), 'LineWidth', 4)
% legend('$|z\widehat{\gamma}(z)|$', 'm', 'M','Interpreter','latex')
xlabel('log_{10}(\tau)')
ylabel('log_{10} of the function value')
tt1 = title('$\alpha_1|\widehat{\gamma}(z)|^2 + \alpha_2|z\widehat{\gamma}(z)|^2$','Interpreter','latex', 'FontSize',20);
% title('The value of $|z\widehat{\gamma}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')
% ylim([])

nexttile
hold on;grid on;
plot(log10(tau_grid_real), log10(h_eps_z.^2), 'LineWidth', 4)
yline(log10(M_h_eps_z), 'LineWidth', 4)
yline(log10(m_h_eps_z), 'LineWidth', 4)
% legend('$|\widehat{h_{\varepsilon}}(z)|^2$', 'm', 'M','Interpreter','latex', 'Location', 'best')
xlabel('log_{10}(\tau)')
tt2 =title('$|\widehat{h_{\varepsilon}}(z)|^2$','Interpreter','latex', 'FontSize',20);
% title('The value of $|\widehat{h_{\varepsilon}}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')
% ylim([-1, y_up])

% subplot(132);
nexttile
hold on;grid on;
plot(log10(tau_grid_real), log10(z_h_eps_z.^2), 'LineWidth', 4)
yline(log10(M_z_h_eps_z), 'LineWidth', 4)
yline(log10(m_z_h_eps_z), 'LineWidth', 4)
% legend('$|z\widehat{h_{\varepsilon}}(z)|^2$', 'm', 'M','Interpreter','latex', 'Location', 'best')
xlabel('log_{10}(\tau)')
tt3 = title('$|z\widehat{h_{\varepsilon}}(z)|^2$','Interpreter','latex', 'FontSize',20);
% title('The value of $|z\widehat{h_{\varepsilon}}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')
% ylim([-1, y_up])



% subplot(133);
nexttile
hold on;grid on;
plot(log10(tau_grid_real), log10(h_eps_z_comb), 'LineWidth', 4)
yline(log10(m_h_eps_z_comb), 'LineWidth', 4)
yline(log10(M_h_eps_z_comb), 'LineWidth', 4)
% legend('$\alpha |z\widehat{h_{\varepsilon}}(z)|^2 + |\widehat{h_{\varepsilon}}(z)|^2$', 'm', 'M','Interpreter','latex', 'Location', 'best')
xlabel('log_{10}(\tau)')
tt4 = title('${\alpha_1 |\widehat{h_{\varepsilon}}(z)|^2 + \alpha_2|z\widehat{h_{\varepsilon}}(z)|^2}$','Interpreter','latex', 'FontSize',20);

% title('The value of $\alpha |z\widehat{h_{\varepsilon}}(z)|^2 + |\widehat{h_{\varepsilon}}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')
% ylim([-1, y_up])


set(findall(gcf,'-property','FontSize'),'FontSize',17)
ttlftsz = 20;
set(tt1, 'fontsize', ttlftsz)
set(tt2, 'fontsize', ttlftsz)
set(tt3, 'fontsize', ttlftsz)
set(tt4, 'fontsize', ttlftsz)


set(fig_coer, 'position', [300, 300, 1000, 300])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

m_h_eps_z_comb
upbd_iprv = (2/m_h_eps_z_comb)*(M_gamma*h_error_paper + g_error_paper);
% sqrt(2*(1/m_h_eps_z_comb*dg_error.^2 + M_gamma_eps_z_comb/m_h_eps_z_comb*h_error.^2));
upbd_iprv

if save_coercivity
    print(fig_coer,  [fig_dir, 'lap_h_gamma.pdf'], '-dpdf', '-r0')
end
