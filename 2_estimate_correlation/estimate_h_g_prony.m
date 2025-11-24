function [EstResult, PronyDetails] = estimate_h_g_prony(I_h, I_f, sysInfo, corr_func, corr_func_obs)
result_h = prony_method(I_h);

F_type  = sysInfo.F_type;
rho     = sysInfo.rho;
dt      = sysInfo.dt;
tgrid   = sysInfo.tgrid;
m       = sysInfo.m;

if sysInfo.mu == 0 &&  (strcmp(F_type, '0') || strcmp(F_type, 'linear'))

    result_g.g = @(x) -m*result_h.dh(x);
    result_g.dg = @(x) -m*result_h.ddh(x);
    result_f = @(x) 0*x;
else
    result_f = prony_method(I_f);
    result_g.g = @(x) -m*result_h.dh(x) + result_f.h(x);
    result_g.dg = @(x) -m*result_h.ddh(x) + result_f.dh(x);
end

PronyDetails.result_h   = result_h;
PronyDetails.result_g   = result_g;
PronyDetails.result_phi = result_f;

EstResult.h    = result_h.h;
EstResult.g    = result_g.g;
EstResult.dg   = result_g.dg;


%%
a2 = real(-sum(result_h.w./result_h.lam));
a1 = result_h.h(0);
comb_coef = a1^2/(a1^2 + a2^2);
% comb_coef is alpha_1 in the paper

%
h_error = sqrt(sum((corr_func.h - EstResult.h(tgrid)).^2.*rho(tgrid))*dt);
g_error = sqrt((1-comb_coef)*sum((corr_func.g - EstResult.g(tgrid)).^2.*rho(tgrid))*dt + comb_coef*sum((corr_func.dg - EstResult.dg(tgrid)).^2.*rho(tgrid))*dt);

h_norm = sqrt(sum((corr_func.h).^2.*rho(tgrid))*dt);
g_norm = sqrt((1-comb_coef)*sum((corr_func.g).^2.*rho(tgrid))*dt + comb_coef*sum((corr_func.dg).^2.*rho(tgrid))*dt);

h_rel_error = h_error/h_norm;
g_rel_error = g_error/g_norm;

EstResult.comb_coef = comb_coef;
EstResult.error.h = h_error;
EstResult.error.g = g_error;
EstResult.error.h_norm = h_norm;
EstResult.error.g_norm = g_norm;
EstResult.error.h_rel = h_rel_error;
EstResult.error.g_rel = g_rel_error;


fprintf('--- Estimation Results ---\n');
fprintf('Combined Coefficient:   %.4f\n', EstResult.comb_coef);
fprintf('h Error (abs):          %.4e\n', EstResult.error.h);
fprintf('g Error (abs):          %.4e\n', EstResult.error.g);
fprintf('h relative Error (abs): %.4e\n', EstResult.error.h_rel);
fprintf('g relative Error (abs): %.4e\n', EstResult.error.g_rel);
%%

Prony_g     = EstResult.g;
Prony_dg    = EstResult.dg;
Prony_h     = EstResult.h;
VACF_obs_len= I_h.prony_N;
tgrid_obs   = I_h.tgrid_obs;

figure;
t = tiledlayout(1, 3, "TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';

scale = max(abs(corr_func.h))*0.2;

nexttile;hold on;grid on;

% fill([tgrid_obs, fliplr(tgrid_obs)], ...
%      [scale*rho(tgrid_obs), zeros(size(tgrid_obs))], ...
%      [0.9, 0.9, 0.9], ...         % light gray fill
%      'EdgeColor', 'none');       % no border
plot(tgrid, corr_func.h, 'LineWidth', 3, 'DisplayName', 'True h')
plot(tgrid, Prony_h(tgrid), '-.', 'LineWidth', 4, 'DisplayName', 'Prony est h')
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.h(1:VACF_obs_len), '-o', 'LineWidth', 2, 'DisplayName', 'Observe h')
if isfield(sysInfo, 'sym')
    plot(tgrid, sysInfo.sym.h_true(tgrid), 'LineWidth', 3, 'DisplayName', 'Symbolic h');
end
xlim([0, min(20, tgrid(end))])

title("h: Auto covariance of v");
legend()
xlabel('time t')

nexttile;hold on;grid on;
% fill([tgrid_obs, fliplr(tgrid_obs)], ...
%      [rho(tgrid_obs)*scale, zeros(size(tgrid_obs))], ...
%      [0.9, 0.9, 0.9], ...         % light gray fill
%      'EdgeColor', 'none');       % no border
plot(tgrid, corr_func.g, 'LineWidth', 3, 'DisplayName', 'True g')
plot(tgrid, Prony_g(tgrid), '-.', 'LineWidth', 4, 'DisplayName', 'Prony est g')
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.g(1:VACF_obs_len), '-o', 'LineWidth', 2, 'DisplayName', 'Observe g')
if isfield(sysInfo, 'sym')
    plot(tgrid, sysInfo.sym.g_true(tgrid), 'LineWidth', 3, 'DisplayName', 'Symbolic g');
end
xlim([0, min(20, tgrid(end))])

title("g: g = -h'");
legend()
xlabel('time t')


nexttile;hold on;grid on;
% fill([tgrid_obs, fliplr(tgrid_obs)], ...
%      [rho(tgrid_obs)*scale, zeros(size(tgrid_obs))], ...
%      [0.9, 0.9, 0.9], ...         % light gray fill
%      'EdgeColor', 'none');       % no border
plot(tgrid, corr_func.dg, 'LineWidth', 3, 'DisplayName', 'True dg')
plot(tgrid, Prony_dg(tgrid), '-.', 'LineWidth', 4, 'DisplayName', 'Prony est dg')
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.dg(1:VACF_obs_len), '-o', 'LineWidth', 2, 'DisplayName', 'Observe dg')
if isfield(sysInfo, 'sym')
    plot(tgrid, sysInfo.sym.dg_true(tgrid), 'LineWidth', 3, 'DisplayName', 'Symbolic dg');
end
xlim([0, min(20, tgrid(end))])



title("g': g' = -h''");
legend()
xlabel('time t')



set(gcf, 'position', [300, 300, 800, 220])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',14)


% print(fig3, [fig_dir, 'h_g_est.pdf'], '-dpdf', '-r0')

end