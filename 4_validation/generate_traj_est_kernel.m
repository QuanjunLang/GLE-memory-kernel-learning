function [v, f, h] = generate_traj_est_kernel(R, S, F, M, para, random_seed, VACF_gen_len, obs_gap)


loadON = para.loadON;
wu = para.wu;
N = para.N;
F_type = para.F_type;
method = para.method;



dt = pi/wu;
dw = wu/N;
T0 = 2*pi/dw;
tgrid = 0:dt:M*dt-dt;
T = M*dt;

id = ['est_traj_data/', num2str(random_seed), '_', F_type, '_', num2str(para.beta), '_', method];
filename = [id, '_traj.mat'];
noisename = [id, '_nse.mat'];


if exist(filename, 'file') == 2 && loadON
    load(filename, 'v', 'f', 'h');
else
    %     generate_correlated_noise
    if exist(noisename, 'file') == 2 && loadON
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
        force = F(v(i-1)) - memory + f(i-1);
        v(i) = v(i-1) + dt*force;
    end
    
end


h = xcov(v, v, 'unbiased')';
h = h(M:end);
save(filename, 'v', 'f', 'h');


end