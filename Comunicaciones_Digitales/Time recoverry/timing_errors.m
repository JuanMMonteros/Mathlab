%%

clear 
close all
clc;

% Parametros 
BR = 32e9;      % Bd
L = 10000;      % Simulation Length
N = 4;          % Oversampling rate
fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 conseutive samples at Tx output
rolloff = 0.1;  % Pulse shaping rolloff 
ps_taps = 10;   % Pulse shaping taps
M=4;            % M-QAM
k=log2(M);

% PARAMETROS DEL CLOCK
clk_phase = 0; % Son fracciones de T: si clk_phase = 0.5 -> t0 =0.5*T=0.5/BR
clk_ppm =  0;  % fs_tx=fs. Relacion fs_rx = fs*(1+ppm*1e-6)
clk_jitter_amp = 0*0.1; % Amplitud de jitter deterministico
clk_jitter_freq = 50e6; % Frecuencia de jitter deterministico 
clk_jitter_std = 5e-12; % Desviacion estandar de jitter random

%% Two symbols generation (+1,-1)
decsymbs = randi([0 M-1], L,1);
x=qammod(decsymbs, M);

% Upsampling to change sampling rate
xup = upsample(x,N);

% Filter to interpolate the signal
h = rcosine(BR, fs, 'sqrt', rolloff, ps_taps);
yup = filter(h, 1, xup); 

%% Errores de timing
Lyup = length(yup);
line = (0:Lyup-1).';

% Tiempo ideal
time_ideal =  line .* Ts;

% Modelo ppm
fs_real = fs * (1+clk_ppm*1e-6);
Ts_real = 1/fs_real;

% Modelo jitter deterministico
Ts_jitter_det = clk_jitter_amp*1/BR * sin(2*pi*clk_jitter_freq*time_ideal);

% Modelo jitter deterministico
Ts_jitter_rand = clk_jitter_std * randn(Lyup, 1);

% Modelo jitter random
time_real = clk_phase*1/BR + line.*Ts_real + Ts_jitter_det + Ts_jitter_rand;

% Aplico los fenomenos resampleando 
yrs = interp1(time_ideal,yup,time_real, 'spline', 0);

%% RECEPTOR
% MF
yrx = filter(h,1,yup);
yrx_rs = filter(h,1,yrs);

% Downsampling
ybd = downsample(yrx,N);
ybd_rs = downsample(yrx_rs,N);

%% PLOTS
eyediagram(yrx(L*N/2:end-1000),2*N);
title('Diagrama de ojo sin ppm')
eyediagram(yrx_rs(L*N/2:end-1000),2*N);
title('Diagrama de ojo con ppm')

figure
subplot(1,2,1)
plot(real(ybd(1:L/2)),'.');
title('Simbolos recibidos sin errores')

subplot(1,2,2)
plot(real(ybd_rs(1:L/2)),'.');
title('Simbolos recibidos con errores')

figure; 
subplot(1,2,1)
plot(time_ideal,'r','Linewidth',2);
hold on;
plot(time_real, '--b','Linewidth',2);
legend('Tiempo ideal','Tiempo real');
title('Tiempo Real vs Ideal')
subplot(1,2,2)
plot(time_ideal-time_real, '-k','Linewidth',2);
title('Diferencia')
    