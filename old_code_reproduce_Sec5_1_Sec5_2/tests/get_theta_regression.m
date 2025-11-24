function [theta, result_details] = get_theta_regression(I, h_est, dg_est)
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
h_conv_phi_dx_method = I.h_conv_phi_dx_method;


xgrid_regression = linspace(lb, rb, N)';

rho_grid_regression = rho(xgrid_regression);
dx = xgrid_regression(2) - xgrid_regression(1);

hgrid_regression = h_est(xgrid_regression);
dg_grid_regression = dg_est(xgrid_regression);

switch I.basis_type
    case 'spline'
        [basis_mat, knots, basis_num] = spline_basis_new(lb, rb, knot_num, deg, deg, xgrid_regression, free_bdry);
        aa = 1;
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


switch h_conv_phi_dx_method
    case 'direct'
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
        
        h_conv_phi_dx = zeros(N, basis_num);
        for k = 1:basis_num
            h_conv_phi_dx(:, k) = gradient(h_conv_phi(:, k), dx);
        end
        
    case 'basis_derivative'
        h_conv_phi_dx = zeros(N, basis_num);
        basis_fun = basis_mat{1};
        basis_fun_dx = basis_mat{2};
        for k = 1:basis_num
            phi_dx = basis_fun_dx(:, k);
            phi = basis_fun(:, k);
            h_conv_phi_dx(:, k) = Volterra_convolution(hgrid_regression, phi_dx, dx);
            h_conv_phi_dx(:, k) = h_conv_phi_dx(:, k) + phi(1)*hgrid_regression;
        end
        a = 1;
end


A = zeros(basis_num, basis_num);
for i = 1:basis_num
    for j = 1:basis_num
        A(i, j) = sum(h_conv_phi_dx(:, i).*h_conv_phi_dx(:, j).*rho_grid_regression)*dx;
    end
end

b = zeros(basis_num, 1);
for i = 1:basis_num
    b(i) = sum(h_conv_phi_dx(:, i).*dg_grid_regression.*rho_grid_regression)*dx;
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
        
        basis_fourier_transform = spline_fourier_transform(I);
        theta_Fourier_transform_temp = @(y) 0*y;
        for j = 1:length(c)
            theta_Fourier_transform_temp = @(y) theta_Fourier_transform_temp(y) + c(j)*basis_fourier_transform{j}(y);
        end
        
        %%
        theta_Fourier_transform_true = @(y) (theta_Fourier_transform_temp(y) + conj(theta_Fourier_transform_temp(y)))/(2*pi);
        result_details.theta_Fourier_transform = theta_Fourier_transform_true;
        result_details.c = c;
        
        
        
    case 'prony'
        theta = @(x) 0*x;
        for i = 1:basis_num
            theta = @(x) theta(x) + c(i)*basis_func{1, i}(x);
        end
end



end