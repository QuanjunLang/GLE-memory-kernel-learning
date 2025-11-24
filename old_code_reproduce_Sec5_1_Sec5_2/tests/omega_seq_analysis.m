

% Using measure for Prony does not make changes
% Use the same for every omega instead
%% Omega sequence
Len = 100;
omega_seq = 10.^linspace(-5, 2, Len);
m_h_eps_seq = zeros(Len, 1);
M_h_eps_seq = zeros(Len, 1);
m_gamma_seq = zeros(Len, 1);
M_gamma_seq = zeros(Len, 1);

theta_regression_error_seq = zeros(Len, 1);
theta_prony_error_seq = zeros(Len, 1);


h_error_seq = zeros(Len, 1);
dg_error_seq = zeros(Len, 1);
gamma_regression_error_seq = zeros(Len, 1);
gamma_prony_error_seq = zeros(Len, 1);

gamma_regression_unif_error_seq = zeros(Len, 1);
gamma_prony_unif_error_seq = zeros(Len, 1);

gamma_regression_ill_error_seq = zeros(Len, 1);
gamma_regression_ill_unif_error_seq = zeros(Len, 1);


prony_use_rho = 0;
if ~prony_use_rho
    Prony_I.rho = rho;
    Prony_I.weight_method = 'LS_freq';
    result_prony = prony_method(Prony_I);
    
    Prony_g = result_prony.g;
    Prony_dg = result_prony.dg;
    Prony_h = result_prony.h;
    
    result_gamma = get_gamma_from_prony_h(Prony_I, result_prony);
    theta_prony = @(x) real(result_gamma.gamma(x));
end


for i = 1:length(omega_seq)
    i
    omega = omega_seq(i);
    rho = @(x) exp(-omega*x)*omega;        % omega
    
%     if prony_use_rho
%         
%         Prony_I.rho = rho;
%         Prony_I.weight_method = 'LS_freq';
%         result_prony = prony_method(Prony_I);
%         
%         Prony_g = result_prony.g;
%         Prony_dg = result_prony.dg;
%         Prony_h = result_prony.h;
%         result_gamma = get_gamma_from_prony_h(Prony_I, result_prony);
%         theta_prony = result_gamma.gamma;
%     end
    h_error_seq(i) = sqrt(sum((RR - Prony_h(tgrid')).^2.*rho(tgrid'))*dt);
    dg_error_seq(i) = sqrt(sum((dg_true - Prony_dg(tgrid')).^2.*rho(tgrid'))*dt);
    
    [m_h_eps_seq(i), M_h_eps_seq(i)] = get_m_M(omega_seq(i), result_prony.h_lap);
    [m_gamma_seq(i), M_gamma_seq(i)] = get_m_M(omega_seq(i), R_lap);
    
    
    h_est = @(x) real(result_prony.h(x));
    g_est = @(x) real(result_prony.g(x));
    dg_est = @(x) real(result_prony.dg(x));
    
    regression_para.omega = omega;
    regression_para.reg_method = 'RKHS';
    theta = get_theta_regression(regression_para, h_est, dg_est);
    
    gamma_regression_error_seq(i) = sqrt(sum((theta(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
    gamma_prony_error_seq(i) = sqrt(sum((theta_prony(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
    
    gamma_regression_unif_error_seq(i) = sqrt(sum((theta(tgrid') - R(tgrid')).^2)*dt);
    gamma_prony_unif_error_seq(i) = sqrt(sum((theta_prony(tgrid') - R(tgrid')).^2)*dt);
    
    theta_ill = get_theta_regression_ill_posed(regression_para, h_est, g_est);
    gamma_regression_ill_error_seq(i) = sqrt(sum((theta_ill(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
    gamma_regression_ill_unif_error_seq(i) = sqrt(sum((theta_ill(tgrid') - R(tgrid')).^2)*dt);
    
end

%%
regression_para.omega = 0;
regression_para.reg_method = 'RKHS';
theta = get_theta_regression(regression_para, h_est, dg_est);
a =  sqrt(sum((theta(tgrid') - R(tgrid')).^2)*dt);
%%
figure; hold on;
plot(log10(omega_seq), log10(m_h_eps_seq), 'o', 'LineWidth', 4);
plot(log10(omega_seq), log10(M_h_eps_seq), 'o', 'LineWidth', 4);
plot(log10(omega_seq), log10(m_gamma_seq), 'LineWidth', 4);
plot(log10(omega_seq), log10(M_gamma_seq), 'LineWidth', 4);
legend('m_h', 'M_h', 'm_\gamma' ,'M_\gamma')

%%
h_norm_omega = zeros(Len, 1);
dg_norm_omega = zeros(Len, 1);
theta_norm_omega = zeros(Len, 1);

for i = 1:Len
    i
    omega = omega_seq(i);
    rho = @(x) exp(-omega*x)*omega;        % omega
    
    h_norm_omega(i) = sqrt(sum(h_true(tgrid).^2.*rho(tgrid)));
    dg_norm_omega(i) = sqrt(sum(dg_true.^2.*rho(tgrid')));
    theta_norm_omega(i) = sqrt(sum(R(tgrid).^2.*rho(tgrid)));
end
h_norm_unif = sqrt(sum(h_true(tgrid).^2));
dg_norm_unif = sqrt(sum(dg_true.^2));
theta_norm_unif = sqrt(sum(R(tgrid).^2));

%%
figure;hold on;
plot(log10(omega_seq), log10(h_error_seq))
plot(log10(omega_seq), log10(dg_error_seq))
%%
figure; hold on; grid on;
error_upper_bound_seq = 1./m_h_eps_seq.*(dg_error_seq + M_gamma_seq.*h_error_seq);
% loglog(omega_seq, upper_bound_errror, 'LineWidth', 4);
plot(log10(omega_seq), log10(error_upper_bound_seq), '-', 'Color', 'green','DisplayName', 'error_upper_bound', 'LineWidth', 4);
plot(log10(omega_seq), log10(1./m_h_eps_seq), '.', 'Color', 'green', 'DisplayName', 'coercivity_constant', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_error_seq), '-', 'Color', 'black', 'DisplayName', 'regression error L^2(rho)', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_prony_error_seq), '-',  'Color', 'blue', 'DisplayName', 'Prony error L^2(rho)', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_unif_error_seq),'-.', 'Color', 'black', 'DisplayName', 'regression error L^2(Lebesgue)', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_prony_unif_error_seq), '-.',  'Color', 'blue', 'DisplayName', 'Prony error L^2(Lebesgue)', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_ill_error_seq),  'Color', 'red', 'DisplayName', 'ill-posed regression error L^2(rho)', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_ill_unif_error_seq), '-.', 'Color', 'red', 'DisplayName', 'ill-posed regression error L^2(Lebesgue)', 'LineWidth', 4);
yline(log10(a))
xlabel('log 10 omega')
ylabel('log 10 L2 error')
legend()
%%

figure; hold on; grid on;
error_upper_bound_seq = 1./m_h_eps_seq.*(dg_error_seq + M_gamma_seq.*h_error_seq);
% loglog(omega_seq, upper_bound_errror, 'LineWidth', 4);
plot(log10(omega_seq), log10(error_upper_bound_seq./theta_norm_omega), '-', 'Color', 'green','DisplayName', 'error_upper_bound', 'LineWidth', 4);
plot(log10(omega_seq), log10(1./m_h_eps_seq./theta_norm_omega), '.', 'Color', 'green', 'DisplayName', 'coercivity_constant', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_error_seq./theta_norm_omega), '-', 'Color', 'black', 'DisplayName', 'regression error L^2(rho)', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_prony_error_seq./theta_norm_omega), '-',  'Color', 'blue', 'DisplayName', 'Prony error L^2(rho)', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_unif_error_seq./theta_norm_unif),'-.', 'Color', 'black', 'DisplayName', 'regression error L^2(Lebesgue)', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_prony_unif_error_seq./theta_norm_unif), '-.',  'Color', 'blue', 'DisplayName', 'Prony error L^2(Lebesgue)', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_ill_error_seq./theta_norm_omega),  'Color', 'red', 'DisplayName', 'ill-posed regression error L^2(rho)', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_ill_unif_error_seq./theta_norm_unif), '-.', 'Color', 'red', 'DisplayName', 'ill-posed regression error L^2(Lebesgue)', 'LineWidth', 4);
yline(log10(a))
xlabel('log 10 omega')
ylabel('log 10 L2 error')
legend()


%%
figure;hold on;
plot(log10(omega_seq), log10(h_error_seq./h_norm_omega))
plot(log10(omega_seq), log10(dg_error_seq./dg_norm_omega))
legend('h rel err', 'dg rel err')
