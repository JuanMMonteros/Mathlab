%%

clear 
close all
clc;

% Parametros 
BR = 32e9;      %Bd
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
clk_ppm =  100;

%% Two symbols generation 
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

% Modelo jitter random
time_real = line.*Ts_real; 

% Aplico los fenomenos resampleando 
yrs = interp1(time_ideal,yup,time_real, 'spline', 0);

%% RECEPTOR
% MF
yrx = filter(h,1,yup);
yrx_rs = filter(h,1,yrs);

% Downsampling
ybd = downsample(yrx,N);
ybd_rs = downsample(yrx_rs,N);
    
%% CORRECCION DE ERRORES DE TIMING

Ntaps = 63;                 % taps del filtro interpolador en simbolos
buffer = zeros(Ntaps,1);
Tint = Ts;                  % Implica TED ideal

interp_out = zeros(L*N,1);  % Salida del interpolador

uk=0;                       % mu
uk_d=0;                     % mu anterior
signo_d = 0;                % signo anterior

uk_log = zeros(L*N,1);      % Para despues graficar uk

extra_mem=0;                % Implica descartar o repetir una muestra
init_timer = Ntaps;

for m = 1:L*N-200
    
    index = m+extra_mem;
    buffer = yrx_rs(index:index+Ntaps);
    
    % Filtro interpolador
    taps = interp_filter(Ntaps, uk); %uk entre 0 y 1
    
    %figure(1);
    %stem(taps); hold on;
    
    % Guardo uk para luego plotear
    uk_log(m)=uk;%*Ts_real;
    
    % Aplico correccion
    yf = conv(taps, buffer);
    yf = yf(Ntaps+1); 
    interp_out(m) = yf;
    
    % Calculo el proximo mu a corregir
    uk = mod(uk + Ts/Ts_real, 1);
    
    % Manejo del desbordamiento de la FIFO (Parte Entera de mu)
    signo = sign(uk-uk_d);
    if m>init_timer
        if signo ~= signo_d % hubo un overflow o underflow
            if signo==1 && signo_d ==-1 % UNDERFLOW
              extra_mem =extra_mem - 1;
              signo=-1; % IMPORTANTE
            else % OVERFLOW
              extra_mem = extra_mem + 1;  
              signo=1; % IMPORTANTE
            end
        end
    end
    
    % Guardo el mu actual 
    uk_d = uk;
    signo_d = signo;
end

%% PLOTS

eyediagram(yrx(L*N/2:end-1000),2*N);
title('Diagrama de ojo sin ppm')
eyediagram(yrx_rs(L*N/2:end-1000),2*N);
title('Diagrama de ojo con ppm')

figure
subplot(1,2,1)
plot(real(ybd(1:L/2)),'.');
title('Simbolos recibidos sin ppm')

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

ybd2 = downsample(interp_out,N,0);
figure
s2plot=ybd2; Lsp = length(s2plot);
t=(0:Lsp-1).*Tint;
plot(real(s2plot),'.');
title('Simbolos recuperados sin ppm')

figure
plot((uk_log(1:end-1000)))
title('Entrada del interpolador');
