%%
clear all
close all

% Primero vamos a armar un ecualizador 
% solo con modo DD (dirigido por decisiones)

% Parametros 
BR = 32e9; %Bd
L = 250000; % Simulati on Length
N = 2; % Oversampling rate
fs = N*BR; % Sampling rate to emulate analog domain
T = 1/BR; % Time interval between two consecutive symbols
Ts = 1/fs; % Time between 2 conseutive samples at Tx output
EbNo = 50; %dB
M=64; %M-QAM
k=log2(M);

% Two symbols generation (+1,-1)
decsymbs = randi([0 M-1], L,1);
x=qammod(decsymbs, M);

% Upsampling to change sampling rate
xup = upsample(x,N);

% Filter to interpolate the signal
h = rcosine(BR,fs,'sqrt', 0.1,10);
yup = filter(conv(h, 1),1,xup); %fir1(9,0.8*1/N)


%% Agregamos ISI y ruido
% 1) Agrego ISI
ch = fir1(8, 0.35); 
ych = filter(ch,1,yup);

% 2) Del EbNo, obtengo la SNR
SNR_dB = EbNo - 10*log10(N) + 10*log10(k);
SNR = 10^(SNR_dB/10);
Pt = var(ych);
No = Pt/SNR;
sigma = sqrt(No/2);
ruido = sigma.*randn(length(ych),1)+1j.*sigma.*randn(length(ych),1);

rx1 = ych+ruido;

h = 1;%fir1(3, 0.45); 
rx = filter(h,1,rx1);

%% AGC
% Mantiene la potencia de entrada del FFE constante
% para evitar que se mueva el paso optimo
target = 0.3; %VRMS
metric = std(rx);
agc_gain = target/metric;
rx_norm = rx.*agc_gain;

%% ECUALIZADOR FSE 
% Vamos a empezar solo con DD (slicer) sin FCR (carrier recovery)
NTAPS = 63; % Por comodidad impar
step = .2e-2;%2e-3;

i_equalizer = rx_norm(1:end); % Columna
Lrx = length(i_equalizer);
buffer = zeros(NTAPS,1); %Columna

o_filter = zeros(Lrx,1); % salida del FSE (tasa 2)
yk = zeros(fix(Lrx/N),1); % salida del FSE (tasa 1)
ak = zeros(fix(Lrx/N),1); % salida del slicer (tasa 1)
W = zeros(NTAPS,1); W((NTAPS+1)/2)=1.0;
tap_leak_gain = 0*1e-2;

for m=1:Lrx
    % Actualizo buffer de entrada
    buffer(2:end)=buffer(1:end-1);
    buffer(1) = i_equalizer(m);
    
    % Filtro
    yfilt = W.'*buffer; %OJO: ' transpone y conjuga! usar .'
    o_filter(m) = yfilt;
    
    if mod(m,N)==0 %Decimo por N
        % Este slicer es lentoooo, cambiarlo por uno mas eficiente
        yk(m/N) = yfilt;
        slicer_out = slicer(yfilt,M);
        ak(m/N)=slicer_out;
        
        % Calculo del error
        Ek = slicer_out-yfilt;
        
        % Emulo el sobre-muestreo del error, ejecutando el LMS
        % dentro de este if
        W = W*(1-step*tap_leak_gain) + step.*Ek.*conj(buffer);
        %W = W + step.*Ek.*conj(buffer);
    end
end


%% Plots de interes
% Respuesta al impulso del FFE
figure
subplot 211
stem(real(W))
title('Parte Real')
subplot 212
stem(imag(W))
title('Parte Imaginaria')

% Respuesta en fcia del FFE
NFFT = 256;
fstep = fs/NFFT;
f = -fs/2:fstep:fs/2-fstep;
figure
plot(f, abs(fftshift(fft(W,NFFT))));
title('Modulo de la rta. en fcia. del FFE')

%%
% Constelacion a la salida del FFE (tasa 1)
% Siempre que hagan algo adaptivo evitar las primeras muestras
scatterplot(yk(L/2:end-100));


%% ALINEO
T0 = L/2;
delay = finddelay(x,ak);
x2 = ak(T0+delay:end);
x1 = x(T0:T0+length(x2)-1);

% SLICER
% Calcula SER simulado (no voy a calcular BER real por comodidad)
% Uds. calculen el BER como debe ser
SER= sum(x2~=x1)/length(x1);
ber_aprox = 1/k*SER
ber_teorico = berawgn(EbNo, 'qam', M)

%%
figure
NFFT=1024;
fv=(0:NFFT-1)/NFFT*fs/1e9;
% Rta canal
H_CANAL = abs(fft(ch,NFFT)).^2;
plot(fv(1:NFFT/2), H_CANAL(1:NFFT/2)/H_CANAL(1))

hold all
% PSD Tx
[P0,f0]=pwelch(yup,hanning(NFFT/4),0,NFFT,fs/1e9);
plot(f0(1:NFFT/2),P0(1:NFFT/2)/mean(P0(1:10)));
% PSD canal
[P0,f0]=pwelch(ych,hanning(NFFT/4),0,NFFT,fs/1e9);
plot(f0(1:NFFT/2),P0(1:NFFT/2)/mean(P0(1:10)));

% Rta FFE
H_FFE = abs(fft(W,NFFT)).^2;
plot(fv(1:NFFT/2), H_FFE(1:NFFT/2)/H_FFE(1))

% PSD salida FFE
[P0,f0]=pwelch(o_filter(end-L:end),hanning(NFFT/4),0,NFFT,fs/1e9);
plot(f0(1:NFFT/2),P0(1:NFFT/2)/mean(P0(1:10)));

grid on
xlabel('Frequency [GHz]')
ylabel('Amplitude')
legend('Rta. Fcia. Canal', 'PSD Tx', 'PSD Canal', 'Rta FFE', 'PSD Salida FFE')


