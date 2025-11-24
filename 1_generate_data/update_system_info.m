function sysInfo = update_system_info(sysInfo, memoryInfo)



%% Assigning the external force term F
if ~strcmp(sysInfo.F_type, 'linear')
    sysInfo.mu = 0;
end


switch sysInfo.F_type
    case '0'
        F = @(x) 0*x;
    case 'linear'
        F = @(x) -sysInfo.mu*x;
    case 'double_well'
        E0 = 0.02;
        F = @(x) -E0*x.*(x.^2 - 1);
    case 'single_well'
        E0 = 0.5;
        F = @(x) -E0*x;
    case 'steep_double_well'
        E0 = 0.1;
        F = @(x) -E0*x.*(x.^2 - 4);
    case 'very_steep_double_well'
        E0 = 1;
        F = @(x) -E0*x.*(x.^2 - 4);
    case 'very_steep_and_far_double_well'
        E0 = 1;
        F = @(x) -E0*x.*(x.^2 - 16);
end

%% Assigning the drift term G

switch sysInfo.G_type
    case '0'
        G = @(x) 0*x;
    case 'Duffing'
        G_gamma = 0.1;
        G_omega = 1;
        G = @(x) G_gamma*cos(G_omega*x);
    case 'Strong_Duffing'
        G_gamma = 1;
        G_omega = 1;
        G = @(x) G_gamma*cos(G_omega*x);
    case 'Extra_Strong_Duffing'
        G_gamma = 10;
        G_omega = 10;
        G = @(x) G_gamma*cos(G_omega*x);
    case 'Extra_Strong_slow_Duffing'
        G_gamma = 10;
        G_omega = 1;
        G = @(x) G_gamma*cos(G_omega*x);
end



sysInfo.F = F;
sysInfo.G = G;


%% Generate the true autocorrelation when F is linear or 0

switch sysInfo.F_type
    case {'0', 'linear'}
        m   = sysInfo.m;
        mu  = sysInfo.mu;
        h0  = memoryInfo.kernel(0);
        if isfield(memoryInfo, 'kernel_Laplace')
            kernel_Laplace = memoryInfo.kernel_Laplace;

            sysInfo.h_Laplace  = @(x) m*h0./(m*x + mu + kernel_Laplace(x));
        end
end

switch sysInfo.F_type
    case {'0', 'linear'}
        if isfield(memoryInfo, 'w_seq')
            mu = sysInfo.mu;
            m = sysInfo.m;
            beta = sysInfo.beta;
            h0 = memoryInfo.kernel(0);h0 = 1;

            w_seq       = memoryInfo.w_seq;
            lambda_seq  = memoryInfo.lambda_seq;


            syms s
            k = length(w_seq);

            R_lap_sym = 0*s;
            for i = 1:k
                R_lap_sym = R_lap_sym + w_seq(i)/(s - lambda_seq(i));
            end

            h_lap_sym = vpa(1/(m*s + mu + R_lap_sym)/beta);
            h_lap_sym = partfrac(h_lap_sym, 'FactorMode', 'complex');
            % pretty(h_lap_sym)
            C = children(h_lap_sym);

            for i = 1:length(C)
                [N, D] = numden(C{i});
                coef = coeffs(D);
                par = coef(2);
                h_w_seq(i) = N/par;
                h_lam_seq(i) = coef(1)/par;
            end
            h_lam_seq = -h_lam_seq;
            % r = exp(h_lam_seq);
            h_w_seq = h_w_seq / sum(h_w_seq);

            h_true = 0*s;
            for i = 1:length(C)
                h_true = h_true + h_w_seq(i)*exp(s*h_lam_seq(i));
            end
            h_true = matlabFunction(h_true);
            h_true = @(x) real(h_true(x));

            dh_true = 0*s;
            for i = 1:length(C)
                dh_true = dh_true + h_lam_seq(i)*h_w_seq(i)*exp(s*h_lam_seq(i));
            end
            dh_true = matlabFunction(dh_true);
            dh_true = @(x) real(dh_true(x));

            ddh_true = 0*s;
            for i = 1:length(C)
                ddh_true = ddh_true + h_lam_seq(i)*h_lam_seq(i)*h_w_seq(i)*exp(s*h_lam_seq(i));
            end
            ddh_true = matlabFunction(ddh_true);
            ddh_true = @(x) real(ddh_true(x));

            %%
            h_true_lap = 0*s;
            for i = 1:length(C)
                h_true_lap = h_true_lap + h_w_seq(i)/(s - h_lam_seq(i));
            end

            h_true_lap = matlabFunction(h_true_lap);
            h_true_lap = @(x) real(h_true_lap(x));

            sysInfo.sym.h_true = h_true;
            sysInfo.sym.g_true = @(x) -m*dh_true(x) - mu * h_true(x);
            sysInfo.sym.dg_true = @(x) -m*ddh_true(x)- mu * dh_true(x);


        end
end


end