function [corr_func, corr_func_obs] = get_all_correlations_all_method(all_v, all_v_obs, dt, dt_obs, F, corr_method, R, tgrid)
%% correlation functions
switch corr_method
    case 'single_traj'
        % 2. using Temporal integrals for the first trajectory
        corr_func_1traj = get_all_correlations_uncentered(all_v(:, 1)', dt, F);
        corr_func = corr_func_1traj;
        corr_func_obs = get_all_correlations_uncentered(all_v_obs(:, 1)', dt_obs, F);
    case 'monte_carlo'
        % 1. using Spacial domain Monte Carlo integral
        time_lag = 5*dt;   % This is the time lag used here. 0 is using purely Spaciel.
                            % Enlarge it gives extra smoothing for the covariance functions
        corr_func_space = get_all_correlations_ensemble(all_v, dt, F, time_lag);
        corr_func = corr_func_space;
        corr_func_obs = get_all_correlations_ensemble(all_v_obs, dt_obs, F, time_lag);
    case 'mixed'
        % 3. using all Temporal and Spacial integrals
        corr_func_mixed = get_all_correlations_ensemble_mixed(all_v, dt, F);
        corr_func = corr_func_mixed;
        corr_func_obs = get_all_correlations_ensemble_mixed(all_v_obs, dt_obs, F);
end

phi_tgrid = R(tgrid)';
corr_func.h_conv_phi = Volterra_convolution(corr_func.h, phi_tgrid, dt);
plot_h_conv_phi(corr_func, tgrid', corr_method)

%% Compare the performance of the covariance functions

% 
% 
% corr_func_space.h_conv_phi = Volterra_convolution(corr_func_space.h, phi_tgrid, dt);
% corr_func_1traj.h_conv_phi = Volterra_convolution(corr_func_1traj.h, phi_tgrid, dt);
% % corr_func_mixed.h_conv_phi = Volterra_convolution(corr_func_mixed.h, phi_tgrid, dt);
% 
% 
% %%
% plot_h_conv_phi(corr_func, tgrid', 'space integral')
% plot_h_conv_phi(corr_func_1traj, tgrid', '1traj integral')
% % plot_h_conv_phi(corr_func_mixed, tgrid', 'mixed integral')
% %%
% % figure;plot(tgrid, sum(all_v.^2, 2)/M)
% 
% %%
% % figure; hold on;
% % plot(tgrid, corr_func_space.df)
% % plot(tgrid, corr_func_space.dh)
% % plot(tgrid, gradient(corr_func_space.dh, dt))
% % plot(tgrid, corr_func_space.ddh)
% 
% %% Choose a correlation method
% 
% 


end