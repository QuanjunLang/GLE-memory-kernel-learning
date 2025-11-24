x = linspace(0, 10, 1000);
f = @(x) pi./(cosh(x*pi/2));

plot(x, f(x));
plot(x, 1./f(x));