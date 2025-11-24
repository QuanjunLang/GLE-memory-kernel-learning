function [all_v, dt, tgrid, T0, all_f] = generate_trajectories_ensemble(R, S, F, M, para, random_seed)

loadON      = para.loadON;
wu          = para.wu;
N           = para.N;
F_type      = para.F_type;
traj_num    = para.traj_num;
beta        = para.beta;
G           = para.G;
G_type      = para.G_type;
saveON      = para.saveON;

dt      = pi/wu;
dw      = wu/N;
T0      = 2*pi/dw;
tgrid   = 0:dt:M*dt-dt;

filename = ['traj_data_ensemble/', num2str(traj_num), '_', num2str(random_seed), '_F_type', F_type, '_beta', num2str(para.beta), 'G_type', G_type, '_traj.mat'];

if traj_num > 1
    noisename = ['traj_data_ensemble/', num2str(traj_num), '_', num2str(random_seed), '_nse.mat'];
    noise_var_name = 'all_f';
elseif traj_num == 1
     noisename = ['traj_data/', num2str(random_seed), '_nse.mat'];
     noise_var_name = 'f';
end


if exist(filename, 'file') == 2 && loadON
    load(filename, 'all_v', 'all_f');
else
    %     generate_correlated_noise
    if exist(noisename, 'file') == 2
        load(noisename, noise_var_name)
        if traj_num == 1
            all_f = f;
        end

    else
        all_f = zeros(M, traj_num);
        for k = 1:traj_num
            k
            all_f(:, k) = generate_correlated_noise(R, S, wu, N, M, dt);
        end
        save(noisename, 'all_f')
    end

    R_tgrid = R(tgrid);
    G_tgrid = G(tgrid);
    % generate the trajectory v
    all_v = zeros(traj_num, M);
    parfor k = 1:traj_num
        k
        f = all_f(:, k) / beta;
        v = zeros(1, M);
        v(1) = randn();
        for i = 2:M
            gamma = R_tgrid(1:i-1);
            v_history = v(1:i-1);
            memory = sum(gamma.*fliplr(v_history))*dt;
            force = F(v(i-1)) - memory + f(i-1) + G_tgrid(i-1);
            v(i) = v(i-1) + dt*force;
        end
        all_v(k, :) = v;
    end
    all_v = all_v';
    if saveON
        save(filename, 'all_v', 'all_f');
    end
    % else
    %     load('traj.mat', 'v')


end




end