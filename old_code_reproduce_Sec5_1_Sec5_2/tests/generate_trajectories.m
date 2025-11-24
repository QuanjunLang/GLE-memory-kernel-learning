function [v, dt, tgrid, T0, f] = generate_trajectories(R, S, F, M, para, random_seed)

loadON = para.loadON;
wu = para.wu;
N = para.N;
F_type = para.F_type;



dt = pi/wu;
dw = wu/N;
T0 = 2*pi/dw;
tgrid = 0:dt:M*dt-dt;

filename = ['traj_data/', num2str(random_seed), '_', F_type, '_', num2str(para.beta), '_traj.mat'];
noisename = ['traj_data/', num2str(random_seed), '_nse.mat'];

% if isfield(para, 'G')
%     G_tgrid = para.G(tgrid);
% else
%     G_tgrid = 0*tgrid;
% end

if exist(filename, 'file') == 2 && loadON
    load(filename, 'v', 'f');
else
%     generate_correlated_noise
    if exist(noisename, 'file') == 2
        load(noisename, 'f')
    else
        f = generate_correlated_noise(R, S, wu, N, M, dt);
        save(noisename, 'f')
    end
    
    f = f / para.beta;
    % generate the trajectory v
    v = zeros(1, M);
    v(1) = randn();
    cutoff = M;
    for i = 2:M
        if i > cutoff
            gamma = R(tgrid(1:cutoff));
            v_history = v(i-cutoff:i-1);
        else
            gamma = R(tgrid(1:i-1));
            v_history = v(1:i-1);
        end
        memory = sum(gamma.*fliplr(v_history))*dt;
        force = F(v(i-1)) - memory + f(i-1) + G_tgrid(i-1);
        v(i) = v(i-1) + dt*force;
    end
    
    save(filename, 'v', 'f');
    % else
    %     load('traj.mat', 'v')
    
    
end




end