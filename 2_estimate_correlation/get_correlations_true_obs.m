function [corr_func, corr_func_obs] = get_correlations_true_obs(sysInfo, trajInfo, obsInfo, temporal_spacial_method)


dt          = trajInfo.dt;
F           = sysInfo.F;
dt_obs      = obsInfo.dt_obs;
rho         = sysInfo.rho;
RANGE       = trajInfo.T0/10;
m           = sysInfo.m;


DROP_OUT        = 0; % drop the first few time steps
DROP_OUT_obs    = DROP_OUT/obsInfo.gap;

tgrid       = sysInfo.tgrid(1:end-DROP_OUT);
tgrid_obs   = obsInfo.tgrid_obs(1:end-DROP_OUT_obs);
v           = trajInfo.v(DROP_OUT+1:end, :);
v_obs       = obsInfo.v_obs(DROP_OUT_obs+1:end, :);


% corr_func           = get_all_correlations(v(:, 1), dt, F);
% corr_func_obs       = get_all_correlations(v_obs(:, 1), dt_obs, F);

corr_func           = get_all_correlations_new(v, dt, F, m, temporal_spacial_method);
corr_func_obs       = get_all_correlations_new(v_obs, dt_obs, F, m, temporal_spacial_method);



%%
if isfield(sysInfo, 'sym')
    % Time steps
    dt      = tgrid(2) - tgrid(1);
    dt_obs  = tgrid_obs(2) - tgrid_obs(1);

    % Evaluate rho on both grids
    rho_vals     = rho(tgrid);
    rho_vals_obs = rho(tgrid_obs);

    % Evaluate true functions on both grids
    h_true_vals     = sysInfo.sym.h_true(tgrid);
    h_true_vals_obs = sysInfo.sym.h_true(tgrid_obs);

    g_true_vals     = sysInfo.sym.g_true(tgrid);
    g_true_vals_obs = sysInfo.sym.g_true(tgrid_obs);

    dg_true_vals     = sysInfo.sym.dg_true(tgrid);
    dg_true_vals_obs = sysInfo.sym.dg_true(tgrid_obs);

    % Compute L2_rho errors
    L2rho_h_est  = sqrt(sum(rho_vals     .* (corr_func.h     - h_true_vals).^2)     * dt);
    L2rho_h_obs  = sqrt(sum(rho_vals_obs .* (corr_func_obs.h - h_true_vals_obs).^2) * dt_obs);

    L2rho_g_est  = sqrt(sum(rho_vals     .* (corr_func.g     - g_true_vals).^2)     * dt);
    L2rho_g_obs  = sqrt(sum(rho_vals_obs .* (corr_func_obs.g - g_true_vals_obs).^2) * dt_obs);

    L2rho_dg_est = sqrt(sum(rho_vals     .* (corr_func.dg    - dg_true_vals).^2)    * dt);
    L2rho_dg_obs = sqrt(sum(rho_vals_obs .* (corr_func_obs.dg - dg_true_vals_obs).^2) * dt_obs);

    % Display
    fprintf("\nL2_rho error (h):   est = %.4f, obs = %.4f\n", L2rho_h_est, L2rho_h_obs);
    fprintf("L2_rho error (g):   est = %.4f, obs = %.4f\n", L2rho_g_est, L2rho_g_obs);
    fprintf("L2_rho error (dg):  est = %.4f, obs = %.4f\n\n", L2rho_dg_est, L2rho_dg_obs);

end


%%
figure;


subplot(131); hold on;
p1 = plot(tgrid, corr_func.h, 'LineWidth', 3);
p2 = plot(tgrid_obs, corr_func_obs.h, '-.', 'LineWidth', 5);
if isfield(sysInfo, 'sym')
    p3 = plot(tgrid, sysInfo.sym.h_true(tgrid), 'LineWidth', 3);
    legend([p1, p2, p3], 'Estimated h', 'Observed h', 'True h', 'Location', 'best');
else
    legend([p1, p2], 'Estimated h', 'Observed h', 'Location', 'best');
end




xlim([0, RANGE])
title('h')

subplot(132); hold on;
p1 = plot(tgrid, corr_func.g, 'LineWidth', 3);
p2 = plot(tgrid_obs, corr_func_obs.g, '-.', 'LineWidth', 5);
if isfield(sysInfo, 'sym')
    p3 = plot(tgrid, sysInfo.sym.g_true(tgrid), 'LineWidth', 3);
    legend([p1, p2, p3], 'Estimated g', 'Observed g', 'True g', 'Location', 'best');
else
    legend([p1, p2], 'Estimated g', 'Observed g', 'Location', 'best');
end


% plot(tgrid, -m*corr_func.dh - sysInfo.mu*corr_func.h, '-.', 'LineWidth', 5);



xlim([0, RANGE])




title('g')

subplot(133); hold on;
p1 = plot(tgrid, corr_func.dg, 'LineWidth', 3);
p2 = plot(tgrid_obs, corr_func_obs.dg, '-.', 'LineWidth', 5);
if isfield(sysInfo, 'sym')
    p3 = plot(tgrid, sysInfo.sym.dg_true(tgrid), 'LineWidth', 3);
    legend([p1, p2, p3], 'Estimated g''', 'Observed g''', 'True g''', 'Location', 'best');
else
    legend([p1, p2], 'Estimated g''', 'Observed g''', 'Location', 'best');
end
xlim([0, RANGE])
title('g''')

sgtitle('Compare the estimated correlation function (discrete)')



end