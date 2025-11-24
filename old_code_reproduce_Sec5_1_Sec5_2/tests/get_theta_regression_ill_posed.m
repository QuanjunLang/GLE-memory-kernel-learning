function theta = get_theta_regression_ill_posed(I, h_est, g_est)
% This is the exploration measure
% omega = 0.05;


lb = I.lb;
rb = I.rb;
knot_num = I.knot_num;
deg = I.deg;
free_bdry = I.free_bdry;
N = I.N;
omega = I.omega;
rho = @(x) exp(-omega*x);        % omega
reg_method = I.reg_method;
% h_conv_phi_dx_method = I.h_conv_phi_dx_method;

xgrid_regression = linspace(lb, rb, N)';

rho_grid_regression = rho(xgrid_regression);
dx = xgrid_regression(2) - xgrid_regression(1);

hgrid_regression = h_est(xgrid_regression);
g_grid_regression = g_est(xgrid_regression);


switch I.basis_type
    case 'spline'
        [basis_mat, knots, basis_num] = spline_basis_new(lb, rb, knot_num, deg, deg, xgrid_regression, free_bdry);
    case 'prony'
        lam_temp = unique(I.lam);
        lam = unique(real(lam_temp) + 1i*abs(imag(lam_temp)));
        basis_num = length(lam);
        basis_mat{1} = zeros(N, basis_num);
        basis_mat{2} = zeros(N, basis_num);
        basis_func = cell(2, basis_num);
        k = 1;
        for i = 1:length(lam)
            z = lam(i);
            if imag(z) == 0
                basis_func{1, k} = @(x) exp(x*z);
                basis_func{2, k} = @(x) z*exp(x*z);
                k = k + 1;
            else
                a = real(z);
                b = imag(z);
                basis_func{1, k} = @(x) exp(x*a).*cos(a*x);
                basis_func{1, k+1} = @(x) exp(x*a).*sin(a*x);
                
                basis_func{2, k} = @(x) exp(x*a).*(a*cos(b*x) - b*sin(b*x));
                basis_func{2, k+1} = @(x) exp(x*a).*(a*sin(b*x) + b*cos(b*x));
                
                k = k + 2;
            end
            
        end
        
        basis_num = length(basis_func);
        for i = 1:basis_num
            basis_mat{1}(:, i) = basis_func{1, i}(xgrid_regression);
            basis_mat{2}(:, i) = basis_func{2, i}(xgrid_regression);
        end
end
h_conv_phi = zeros(N, basis_num);
basis_fun = basis_mat{1};
for k = 1:basis_num
    phi = basis_fun(:, k);
    conv_temp = zeros(N, 1);
    %     conv_temp(1) = hgrid_regression(1)*phi(1)*dx/2;
    conv_temp(1) = 0;
    for n = 2:N
        val = (hgrid_regression(1)*phi(n) + hgrid_regression(n)*phi(1))/2;
        val = val + sum(hgrid_regression(2:n-1) .* phi(n-1:-1:2));
        conv_temp(n) = val*dx;
    end
    h_conv_phi(:, k) = conv_temp;
end



A = zeros(basis_num, basis_num);
for i = 1:basis_num
    for j = 1:basis_num
        A(i, j) = sum(h_conv_phi(:, i).*h_conv_phi(:, j).*rho_grid_regression)*dx;
    end
end

b = zeros(basis_num, 1);
for i = 1:basis_num
    b(i) = sum(h_conv_phi(:, i).*g_grid_regression.*rho_grid_regression)*dx;
end

switch reg_method
    case {'RKHS', 'ID'}
        [~, c] = L_curve(A, b, 'RKHS', 0);
    case 'LS'
        c = A\b;
        
end

switch I.basis_type
    case 'spline'
sp = spmak(knots, c');
theta = @(x) fnval(sp, x);
    case 'prony'
        theta = @(x) 0*x;
        for i = 1:basis_num
            theta = @(x) theta(x) + c(i)*basis_func{1, i}(x);
        end
end
end