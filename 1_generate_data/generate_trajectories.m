function trajInfo = generate_trajectories(memoryInfo, sysInfo, trajInfo, random_seed)

kernel           = memoryInfo.kernel;
kernel_Fourier   = memoryInfo.kernel_Fourier;
kernel_type      = memoryInfo.kernel_type;


loadON      = trajInfo.loadON;
plotON      = trajInfo.plotON;
wu          = trajInfo.wu;
N           = trajInfo.N;
dt          = trajInfo.dt;
L           = trajInfo.L;
M           = trajInfo.M;

dw          = wu/N;
T0          = 2*pi/dw;
tgrid       = (0:dt:L*dt-dt)';
T           = dt*L;


F_type      = sysInfo.F_type;
G_type      = sysInfo.G_type;
F           = sysInfo.F;
G           = sysInfo.G;
beta        = sysInfo.beta;
m           = sysInfo.m;

data_dir    = sysInfo.data_dir;
fig_dir     = sysInfo.fig_dir;

%% Design the file name system
% Define key-value pairs as a cell array
trj_info_pairs = {
    'KernelType', kernel_type;
    'FType', F_type;
    'GType', G_type;
    'seed', random_seed;
    'beta', beta;
    'm', m;
    'wu', wu;
    'L', L;
    'M', M;
    'dt',dt;
    'N',N;
    'mu', sysInfo.mu;
    };

% Convert to string parts: "key_value"
trj_info_strs = cellfun(@(k,v) [k, '_', num2str(v)], trj_info_pairs(:,1), trj_info_pairs(:,2), 'UniformOutput', false);

% Final filename
trj_file_name = fullfile(data_dir, [strjoin(trj_info_strs, '_'), '_trj.mat']);

%% Noise file name

% Define key-value pairs as a cell array
nse_info_pairs = {
    'KernelType', kernel_type;
    'seed', random_seed;
    'm', m;
    'wu', wu;
    'L', L;
    'M', M;
    'dt',dt;
    'N',N;
    };

% Convert to string parts: "key_value"
nse_info_strs = cellfun(@(k,v) [k, '_', num2str(v)], nse_info_pairs(:,1), nse_info_pairs(:,2), 'UniformOutput', false);

% Final filename
nse_file_name = fullfile(data_dir, [strjoin(nse_info_strs, '_'), '_nse.mat']);


%%
fprintf('Correlated noise fake period: %.2f\n', T0);

if exist(trj_file_name, 'file') == 2 && loadON
    load(trj_file_name, 'f', 'v');
    fprintf('Load trajecory and noise \n')
else

    % generate_correlated_noise
    if exist(nse_file_name, 'file') == 2 && loadON
        load(nse_file_name, 'f')
        fprintf('Load noise \n')
    else
        tic;
        f = zeros(L, M);
        parfor k = 1:M
            k
            f(:, k) = generate_correlated_noise(kernel, kernel_Fourier, wu, N, L, 'dt', dt);
        end
        save(nse_file_name, 'f')
        nse_time = toc;
        fprintf('Generate noise: \t\t%2f second\n', nse_time)
    end

    % generate the trajectory v
    tic;

    f_old = f;
    G_tgrid = G(tgrid);
    f = f / sqrt(beta);

    v = zeros(L, M);
    parfor k = 1:M
        k
        v_temp = zeros(L, 1);
        v_temp(1) = randn();
        cutoff = L;         % Cut off is necessary
        for i = 2:L
            if i > cutoff
                gamma = kernel(tgrid(1:cutoff));
                v_history = v_temp(i-cutoff:i-1);
            else
                gamma = kernel(tgrid(1:i-1));
                v_history = v_temp(1:i-1);
            end
            memory = sum(gamma.*flip(v_history))*dt;
            force = F(v_temp(i-1)) - memory + f(i-1) + G_tgrid(i-1);
            v_temp(i) = v_temp(i-1) + (1/m)*dt*force;
        end
        v(:, k) = v_temp;
    end
    save(trj_file_name, 'v', 'f');
    trj_time = toc;
    fprintf('Generate trajecory: \t\t%2f second\n', trj_time)
end

if plotON
    plot_noise(kernel, kernel_Fourier, tgrid, f(:, 1), wu);
end

trajInfo.f      = f;
trajInfo.v      = v;
trajInfo.T0     = T0;
trajInfo.T      = T;
trajInfo.tgrid  = tgrid;

end