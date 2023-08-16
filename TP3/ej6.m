clc; close all; clear
%%
Lsim=1000e3;
sigma2=0.1;
X = randi([0 1], Lsim,1)*2-1; % Lsim es el n´umero de s´?mbolos a transmitir
N = sqrt(sigma2).*randn(Lsim,1); % sigma2 es la varianza del ruido
Y=X+N; % Se~nal recibida, son los s´?mbolos de TX contaminados con ruido
% Toma de decision en el receptor
Y_decidido = Y;
Y_decidido(Y>0)=1;
Y_decidido(Y<0)=-1;
BER(X~=Y_decidido)=1;
BER = sum(BER);
 