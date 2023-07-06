clear all
close all

sigma1=1;
sigma2=2;

%x1 = sigma.*randn(100e3,1)+mu;
% Genero una gauasiana de varianza sigma^2 y media mu
x1 = sigma1*randn(1000e3,1)+2;
x2 = sigma2*randn(1000e3,1)-5; % var 4

z = x1+x2;

fprintf("Varianza de x1: %2.2f\n", var(x1))
fprintf("Varianza de x2: %2.2f\n", var(x2))
fprintf("Varianza de z: %2.2f\n", var(z))

