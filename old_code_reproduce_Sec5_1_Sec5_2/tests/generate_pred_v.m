function pred_v = generate_pred_v(theta, M, dt, tgrid, F, correlated_noise, initial)

pred_v = zeros(1, M/2);
pred_v(1) = initial;
for i = 2:M/2
    gg = theta(tgrid(1:i-1));
    v_history = pred_v(1:i-1);
    memory = sum(gg.*fliplr(v_history))*dt;
    force = F(pred_v(i-1)) - memory + correlated_noise(i-1);
    pred_v(i) = pred_v(i-1) + dt*force;
end

