regression_para.lb = 0;
regression_para.rb = 30;
regression_para.knot_num = 30;
regression_para.deg = 3;
regression_para.free_bdry = 1;
regression_para.N = 2001;
regression_para.omega = omega;
regression_para.reg_method = 'RKHS';
regression_para.h_conv_phi_dx_method = 'direct';
% regression_para.basis_type = 'prony';
regression_para.basis_type = 'spline';
regression_para.lam = result_gamma.lam_seq;



a1 = sum(result_h.h(tgrid))*dt;
a2 = result_h.h(0);
comb_coef = a1^2/(a1^2 + a2^2);
% comb_coef = 0.5;

regression_para.comb_coef = comb_coef;

h_est = @(x) real(Prony_h(x));
dg_est = @(x) real(Prony_dg(x));
g_est = @(x) real(Prony_g(x));

theta = get_theta_regression(regression_para, h_est, dg_est);
theta_ill_posed = get_theta_regression_ill_posed(regression_para, h_est, g_est);
theta_combined = get_theta_regression_combined(regression_para, h_est, dg_est, g_est);


xgrid_regression = linspace(regression_para.lb, regression_para.rb, regression_para.N)';
dx_regression = xgrid_regression(2) - xgrid_regression(1);

figure;hold on;
plot(xgrid_regression, theta(xgrid_regression), 'LineWidth', 4);
plot(xgrid_regression, theta_ill_posed(xgrid_regression), 'LineWidth', 4);
plot(xgrid_regression, theta_combined(xgrid_regression), 'LineWidth', 4);
plot(xgrid_regression, R(xgrid_regression), 'LineWidth', 4);
legend('Derivative', 'Spacial', 'Combined' ,'True')

gamma_combined_error = sqrt(sum((theta_combined(tgrid_short') - R(tgrid_short')).^2.*rho(tgrid_short'))*dt);
% gamma_combined_error_2 = gamma_combined_error^2

%%
Len = 100;
omega_seq = 10.^linspace(-10, 2, Len);

comb_coef = 0.1;
h = @(x) result_h.h_lap(x);
for i = 1:length(omega_seq)
    i
    omega = omega_seq(i);
    rho = @(x) exp(-omega*x)*omega;        % omega 
    
    tau_grid_real = 10.^linspace(-20, 20, 100000);
    tau_grid = omega + 1i*tau_grid_real;
    hz = abs(h(tau_grid));
    zhz = abs(tau_grid.*h(tau_grid));
    
    m_hz_seq(i) = min(hz)^2;
    M_hz_seq(i) = max(hz)^2;
    
    m_zhz_seq(i) = min(zhz)^2;
    M_zhz_seq(i) = max(zhz)^2;
    
    m_comb_seq(i) = min((1 - comb_coef)*hz.^2 + comb_coef*zhz.^2);
    M_comb_seq(i) = max((1 - comb_coef)*hz.^2 + comb_coef*zhz.^2);
    
    
  
end

figure;hold on;
plot(log10(omega_seq), log10(m_zhz_seq), 'linewidth', 4);
plot(log10(omega_seq), log10(m_hz_seq), 'linewidth', 4);
plot(log10(omega_seq), log10(m_comb_seq), 'linewidth', 4);
legend('zhz', 'hz', 'comb')
xlabel('log_{10}(\omega)')
ylabel('log_{10}(min^2)')


