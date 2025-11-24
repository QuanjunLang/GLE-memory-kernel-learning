function C = compute_corr(A, B, method)
% Computes correlation <A(·), B(·+t)> using method 1 or 2
% A, B: [L x 1] vectors;

L = length(A);

switch method
    case 'xcov'
        C_temp = xcov(A, B, 'unbiased');
        C = C_temp(L:end);
    case 'manual'
        max_lag = L - 1;
        C = zeros(max_lag + 1, 1);

        for tau = 0:max_lag
            total = 0; count = 0;
            if tau < L
                a = A(1:L-tau);
                b = B(1+tau:L);
                total = total + sum(a .* b);
                count = count + (L - tau);
            end
            C(tau+1) = total / count;
        end
    case 'M_only'
        C = zeros(L, 1);
        for tau = 1:L
            C(tau) = mean(A(tau,:) .* B(1,:));  % V(0), V(t)
        end

end
end
