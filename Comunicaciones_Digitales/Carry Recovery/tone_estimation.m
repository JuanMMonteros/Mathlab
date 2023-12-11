%close all
clear all

%% Parametros del simulador
fs = 16e9; 
Ts = 1/fs;
f0 = 200e6;
T0 = 1/f0;
w0 = 2*pi*f0;
NCYCLES = 2000;
tend = NCYCLES*T0;
t = 0:Ts:tend-Ts;

%% Generar sin
A0 = 0.5;
xsin = A0*sin(w0.*t);

fprintf('Valor de pico: %2.2e \n', A0)
fprintf('Valor RMS: %2.2e, STD: %2.2e \n', A0/sqrt(2), std(xsin))
fprintf('Potencia: %2.2e, VAR: %2.2e \n',A0^2/2, var(xsin))

%% Generar ruido 
sigma = sqrt(1);
L = length(xsin);
noise = sigma.*randn(1,L);
fprintf('Noise STD: %2.2e \n', std(noise))
fprintf('Noise VAR: %2.2e \n', var(noise))

%% Sumo senial y ruido
x = xsin+noise;
fprintf('Total VAR: %2.2e \n', var(x))
% Notar que el ruido es aditivo en potencia
figure
plot(x)

%% Calculo del welch
% TAREA: plotear el welch del seno, del ruido y de la suma
NFFT = 10240*16;
%[Pxx, f] = pwelch(x,hanning(NFFT/2),0, NFFT, fs);
[Pxx, f] = pwelch(x,[],0, NFFT, fs);
figure
plot(f,Pxx)

%% Calculo de la potencia a partir del welch
[~, binref] = max(Pxx);
deltaf = fs/NFFT;
req_int_bw = floor(1e6/deltaf)
integral=sum( Pxx(binref-req_int_bw: binref+req_int_bw) );

potencia = integral*deltaf;

A0_est = sqrt(2*potencia);
fprintf('A0 Estimado: %2.2e, A0 Original: %2.2e \n', A0_est, A0) 




