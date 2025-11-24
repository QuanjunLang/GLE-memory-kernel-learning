function obsInfo = generate_observed_trajectory(obsInfo, trajInfo)

% Extract observation parameters
gap   = obsInfo.gap;     % Sampling gap (observation every 'gap' time steps)
std   = obsInfo.std;     % Observation noise standard deviation
len   = obsInfo.len;     % Number of observations

% Extract full trajectory data
tgrid = trajInfo.tgrid;  % Full time grid
v     = trajInfo.v;      % True velocity trajectory

[L, M] = size(v);


% Compute derived observation parameters
dt        = tgrid(2) - tgrid(1);        % Time step size
dt_obs    = gap * dt;                   % Time step for observations
T_obs     = len * dt_obs;               % Total observation time
tgrid_obs = tgrid(1:gap:len*gap);       % Observation time grid
v_obs     = v(1:gap:len*gap, :) + randn(len, M) * std;  % Noisy observed velocities
M_obs     = len;                        % Number of observations

% Store outputs back into obsInfo struct
obsInfo.tgrid_obs = tgrid_obs;
obsInfo.v_obs     = v_obs;
obsInfo.dt_obs    = dt_obs;
obsInfo.T_obs     = T_obs;
obsInfo.M_obs     = M_obs;

% Plotting: show true vs. observed velocity
if obsInfo.plotON
    ind = min(10, M);
    figure('Color', 'w'); hold on;

    % Truncate to observation window
    tgrid_short = tgrid(1:gap*len);
    v_short     = v(1:gap*len, 1:ind);

    % Plot true velocity in light blue
    plot(tgrid_short, v_short, '-', ...
        'Color', [0.5 0.7 1], 'LineWidth', 1.2, ...
        'DisplayName', 'True $v(t)$');

    % Plot sparse noisy observations with smaller markers
    plot(tgrid_obs, v_obs(:, 1:ind), 'o', ...
        'MarkerSize', 4, ...
        'MarkerEdgeColor', [1, 0.4, 0.4], ...
        'MarkerFaceColor', [1, 0.4, 0.4], ...
        'DisplayName', 'Observed $v(t_k)$');

    % Formatting
    xlabel('Time $t$', 'Interpreter', 'latex');
    ylabel('Velocity $v(t)$', 'Interpreter', 'latex');
    title('True and Noisy Observations of Velocity', 'Interpreter', 'latex');

    legend('Location', 'northeast', 'Interpreter', 'latex');
    grid on; box on;

    xlim([0, min(10, tgrid(end))]);
end
end