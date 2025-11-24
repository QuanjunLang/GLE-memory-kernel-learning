function output_seq = evaluation(input_seq, M, Y, phi)


% input_seq = -10:0.1:10;
% input_seq = Y;
output_seq = zeros(M, 1);
for i = 1:length(input_seq)
z = input_seq(i);
result = eye(2, 2);
for j = 1:M
    temp = zeros(2, 2);
    temp(1, 1) = (z - Y(j))/(z - conj(Y(j)));
    temp(1, 2) = phi(j);
    temp(2, 1) = conj(phi(j))*(z - Y(j))/(z - conj(Y(j)));
    temp(2, 2) = 1;
    result = result * temp;
end

param = 0;
theta = (result(1, 1) *param + result(1, 2) ) / (result(2, 1) * param + result(2, 2));
val = 1i*(1+theta)/(1-theta);
output_seq(i) = val;
end