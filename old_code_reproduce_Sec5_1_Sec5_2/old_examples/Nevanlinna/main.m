clc
close all
clear all
rng(1)

% This code generate the Nevanlinna analytic continuation from the paper
% "Nevanlinna Analytical Continuation, Jiani Fei, CHia-Nan Yeh, Emanuel Gull

% We assume a discrete observation of the greens function fermion Matsubara
% points (Matsubara Green's function) and interpolate it along the real
% axis.

%%
% g = @(x) exp(-(x.^2)/2)/sqrt(2*pi) + exp(-((x-5).^2)/2)/sqrt(2*pi);
g = @(x) (exp(-(x.^2)/2)/sqrt(2*pi)+ exp(-((x-1).^2)/2)/sqrt(2*pi));
beta = 100;
M = 36;


Y = 1i*(2*(0:1:M-1)+1)*pi/beta;
C = g(Y);
Lambda = (C-1i)./(C+1i);


%%
% theta = cell(M+1, 1);
% theta{M+1} = @(x) 0*x;
%
% phi = zeros(M, 1);
% phi(1) = C(1);
%

for i = 2:M
    b = com(Y(i), Y(i-1));
    phi(i) = (Lambda(i-1) - Lambda(i))/((Lambda(i)*conj(Lambda(i-1)) - 1)*b);
end


% for i = 1:M
%     j = M-i+1;
%     theta{j} = 1;
% end

%%
phi = zeros(M, 1);
abcds = cell(M, 1);
phi(1) = Lambda(1);



for k = 1:M
    abcds{k} = eye(2, 2);
end


for j = 1:M-1
    for k = j:M
        temp = ones(2, 2);
        temp(1, 1) = (Y(k) - Y(j))/(Y(k) - conj(Y(j)));
        temp(1, 2) = phi(j);
        temp(2, 1) = conj(phi(j))*(Y(k) - Y(j))/(Y(k) - conj(Y(j)));
        temp(2, 2) = 1;
        abcds{k} = abcds{k} * temp;
    end
    phi(j+1) = (-abcds{j+1}(2, 2)*Lambda(j+1) + abcds{j+1}(1, 2))/(abcds{j+1}(2, 1)*Lambda(j+1) - abcds{j+1}(1, 1));
end


%%
input_seq = -5:0.01:5;
output_seq = evaluation(input_seq, M, Y, phi);
figure;hold on;
plot(real(input_seq), real(output_seq));
plot(real(input_seq), g(input_seq));
ylim([min(real(g(input_seq)))-2, max(real(g(input_seq)))+2])

%%
input_seq = Y;
output_seq = evaluation(input_seq, M, Y, phi);
figure;hold on;
plot(imag(input_seq), real(output_seq));
plot(imag(input_seq), g(input_seq));
