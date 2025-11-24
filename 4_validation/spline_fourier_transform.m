function basis_fourier_transform = spline_fourier_transform(I)


% I = regression_para;
lb = I.lb;
rb = I.rb;
knot_num = I.knot_num;
deg = I.deg;
free_bdry = I.free_bdry;
N = I.N;
omega = I.omega;
% rho = @(x) exp(-omega*x);        % omega
% reg_method = I.reg_method;
% h_conv_phi_dx_method = I.h_conv_phi_dx_method;

xgrid_regression = linspace(lb, rb, N)';


[basis_mat, knots, basis_num] = spline_basis_new(lb, rb, knot_num, deg, deg, xgrid_regression, free_bdry);

%%
syms r
% knots = [0, 0, 1];
% deg = 1;
% basis_num = 1;

k = deg+1;

syms t [1 k+1]
syms x y

w = prod(x-t);
wp = diff(w, x);

basis_temp = x-x;
for j = 0:k
    tj = t(j+1);
    basis_temp = basis_temp + exp(-1i*tj*y)/(subs(wp, x, tj));
end
basis_temp = basis_temp*factorial(k)/((-1i*y)^k);



basis_fourier_transform = cell(basis_num, 1);
for j = 1:basis_num
    current_knots = knots(j:j+k);
    if length(unique(current_knots)) == k+1
        temp = subs(basis_temp, t, current_knots);
        basis_fourier_transform{j} = matlabFunction(temp);
    else
        if current_knots(1) == current_knots(2)
            required_len = k+1 - length(unique(current_knots));
            current_knots = r*[(-required_len:1:0)'; zeros(k - required_len, 1)]' + current_knots;
            
        else
            required_len = k+1 - length(unique(current_knots));
            current_knots = r*[zeros(k - required_len, 1); (0:1:required_len)'; ]' + current_knots;
            
        end
%         current_knots
        temp = subs(basis_temp, t, current_knots);
        temp = temp*(current_knots(end) - current_knots(1))/k;
        temp = limit(temp, r, 0);
        basis_fourier_transform{j} = matlabFunction(temp);
    end
end

% pretty(temp)
end

