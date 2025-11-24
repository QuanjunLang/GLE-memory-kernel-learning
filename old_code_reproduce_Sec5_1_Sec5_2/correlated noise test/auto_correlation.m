function RR = auto_correlation(v, dt)

M = length(v);
T = M*dt;
RR = zeros(M, 1);
for t = 1:M
    RR(t) = sum(v(1:end-t+1).*v(t:end))*dt/T;
end

end