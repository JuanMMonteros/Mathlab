clc; close all; clear
%%
Lsim=1e6;%Numero de Simbolos enviados 
sigma2=0.1;% Potencia del ruido
%prubas en las que aumente SNR Po=0.05, dismuni la SNR =0.5
SNR= 1/sigma2;
SNR_db=10*log10(SNR);
X = randi([0 1], Lsim,1)*2-1; % Lsim es el n´umero de s´?mbolos a transmitir
N = sqrt(sigma2).*randn(Lsim,1); % sigma2 es la varianza del ruido
SNR_practica=10*log10(sum(X.*X)/sum(N.*N));
Y=X+N; % Se~nal recibida, son los s´?mbolos de TX contaminados con ruido
% Toma de decision en el receptor
Y_decidido = Y;
Y_decidido(Y>0)=1;
Y_decidido(Y<0)=-1;
F(X~=Y_decidido)=1;%Errores
F = sum(F);
BER= F/Lsim; % bit error rate
BERdb= 10*log10(BER);
ber_teorica = berawgn(SNR_db, 'pam', 2);
ber_teorica = 10*log10(ber_teorica);
fprintf('\n - BER : %.2f db', BERdb)
fprintf('\n - BER Teorica : %.2f db', ber_teorica)
 %% Crear histograma de X
figure;
histogram(X, 'BinMethod', 'auto'); % 'BinMethod' puede ajustarse según tus necesidades
title('Histograma de la señal Digital X');
xlabel('Valores');
ylabel('Frecuencia');

% Crear histograma de N
figure;
histogram(N, 'BinMethod', 'auto');
title('Histograma del Ruido N');
xlabel('Valores');
ylabel('Frecuencia');

% Crear histograma de Y
figure;
histogram(Y, 'BinMethod', 'auto');
title('Histograma de la señal recibida Y');
xlabel('Valores');
ylabel('Frecuencia');