function a = temp_func(s)

f = @(x, s) exp(-s*x) .*(1-3.^x.^2)./ (1 + x.^2);

% s = 1i;  % Example: Laplace evaluated at s = 1

% Use integral with upper bound truncated (e.g., 100)
a = integral(@(x) f(x, s), 0, 100, 'RelTol',1e-10,'AbsTol',1e-12);

end