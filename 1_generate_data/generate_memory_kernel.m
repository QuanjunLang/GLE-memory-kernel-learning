function memoryInfo = generate_memory_kernel(memoryInfo)
% Extract parameters p and q from the structure para

kernel_type = memoryInfo.kernel_type;

switch true
    case contains(kernel_type, 'Test')
        gamma           = @(x) exp(-x);
        gamma_Fourier   = @(x) -1/pi ./ (x.^2 + 1);
        gamma_Laplace   = @(x) 1./(x - 1);

        memoryInfo.w_seq            = 1;
        memoryInfo.lambda_seq       = -1;

    case contains(kernel_type, 'RandomProny')
        % Extract p and q
        p = str2double(regexp(kernel_type, '(?<=_p_)\d+', 'match', 'once'));
        q = str2double(regexp(kernel_type, '(?<=_q_)\d+', 'match', 'once'));
        fprintf('Kernel Type: Random Prony with p = %d complex pairs and q = %d real modes\n', p, q);

        % Initialize sequences for weights and exponents
        w_seq = zeros(2*p+q, 1);
        lambda_seq = zeros(2*p+q, 1);

        % Generate complex conjugate pairs for oscillatory components
        for k = 1:p
            w_seq(2*k-1) = rand();                  % Random weight between -0.5 and 0.5
            lambda_seq(2*k-1) = (rand()-1)*1i - 3*rand(); % Random complex exponent: real part negative
            w_seq(2*k) = w_seq(2*k-1);                      % Conjugate pair has the same weight
            lambda_seq(2*k) = conj(lambda_seq(2*k-1));       % Conjugate exponent
        end

        % Generate purely real, exponentially decaying terms
        for k = 2*p+1:2*p+q
            w_seq(k) = rand();           % Positive random weight
            lambda_seq(k) = -3*rand();   % Negative real exponent
        end

        % Define R(x) and S(x) as zero-initialized anonymous functions
        gamma = @(x) 0*x;
        gamma_Fourier = @(x) 0*x;

        % Build up R(x) and S(x) as sums of exponentials
        for k = 1:2*p+q
            gamma = @(x) gamma(x) + w_seq(k)*exp(lambda_seq(k)*x);           % Sum of weighted exponentials
            gamma_Fourier = @(x) gamma_Fourier(x) + w_seq(k)*(-lambda_seq(k))/pi ./ (x.^2 + lambda_seq(k)^2); % Hilbert transform component
        end

        % Normalize R(x) at x = 0 to be positive
        sgn = sign(gamma(0));
        gamma = @(x) gamma(x)*sgn;
        gamma_Fourier = @(x) gamma_Fourier(x)*sgn;

        % S is the Fourier transform of R
        % The fourier transform used here is
        % F[R] = 1/2pi int R(x) exp(-iwx)dx = S(w)

        gamma_Laplace = @(x) 0*x;
        for k = 1:2*p+q
            gamma_Laplace = @(x) gamma_Laplace(x) + w_seq(k)./(x - lambda_seq(k));
        end
        gamma_Laplace = @(x) gamma_Laplace(x)*sgn;

        memoryInfo.p = p;
        memoryInfo.q = q;
        memoryInfo.numModes         = 2*p+q;
        memoryInfo.kernel_Laplace   = gamma_Laplace;
        memoryInfo.w_seq            = w_seq;
        memoryInfo.lambda_seq       = lambda_seq;

    case contains(kernel_type, 'PowerLaw')
        % Extract a and b
        % tokens = regexp(kernel_type, 'PowerLaw_a_(-?\d+)_b_(-?\d+)', 'tokens');
        % a = str2double(tokens{1}{1});
        % b = str2double(tokens{1}{2});
        fprintf('Kernel Type: PowerLaw (1-3*x.^2)./((1 + x.^2).^3\n');

        gamma = @(x) (1-3*x.^2)./((1 + x.^2).^3);
        gamma_Fourier = @(w) 1/4*w.^2.*exp(-abs(w));
        gamma_Laplace = @(w) -0.5*(w.^2.*(cosint(w) .* sin(w) + (pi/2 - sinint(w)).* cos(w)) - w);
    case contains(kernel_type, 'Polynomial')
        fprintf('Kernel Type: Polynomial (1+x^2)^(-1)\n');

        gamma           = @(x) 1./(1 + x.^2);
        gamma_Fourier   = @(w) 1/2.*exp(-abs(w));
        gamma_Laplace   = @(w) cosint(w) .* sin(w) + (pi/2 - sinint(w)).* cos(w);
end

if memoryInfo.plotON
    range = 30;
    figure;
    subplot(121)
    fplot(gamma, [0, range], 'LineWidth', 4);title('True Kernel \gamma')
    subplot(122)
    fplot(gamma_Fourier, [0, range], 'LineWidth', 4);title('Fourier transform of \gamma')
end



%%
memoryInfo.kernel          = gamma;
memoryInfo.kernel_Fourier  = gamma_Fourier;
memoryInfo.kernel_Laplace   = gamma_Laplace;

end