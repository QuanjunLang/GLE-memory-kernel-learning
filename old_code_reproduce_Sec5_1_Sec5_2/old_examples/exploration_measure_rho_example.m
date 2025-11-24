knots = [1,2,3,4,5,5.5,7,7.2, 9, 10.1, 11, 15];
N = length(knots);

rho = [0];

for i = 1:N-1
    rho = [rho, knots(i+1:end) - knots(i)];
end

rho(1) = [];
rho = rho';

histogram(rho, 100)
