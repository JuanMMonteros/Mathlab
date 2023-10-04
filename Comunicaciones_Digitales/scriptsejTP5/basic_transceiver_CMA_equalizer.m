%%
clear all
close all

% Primero vamos a armar un ecualizador 
% solo con modo DD (dirigido por decisiones)

% Parametros 
BR = 32e9; %Bd
L = .5e6; % Simulati on Length
N = 2; % Oversampling rate RX
Nch = 2; %Oversampling rate TX
fch = Nch*BR; % Sampling rate equalizer
fs = N*BR; % Sampling rate equalizer
T = 1/BR; % Time interval between two consecutive symbols
Ts = 1/fs; % Time between 2 conseutive samples at Tx output
EbNo = 15; %dB
M=16; %M-QAM
k=log2(M);

% Two symbols generation (+1,-1)
decsymbs = randi([0 M-1], L,1);
x=qammod(decsymbs, M);

% RCMA lo calculo aca para que cambie con la modulacion
% automaticamente. Sin embargo, este numero no depende de 
% de los datos, sino de su estadistica. Por ejemplo: para 
% QAM16 uniforme es 3.64
RCMA = sqrt(mean(abs(x(1:1000)).^4)/mean(abs(x(1:1000)).^2));

% Upsampling to change sampling rate
xup = upsample(x,Nch);

% Filter to interpolate the signal
hrc = rcosine(BR,fch,'sqrt', 0.1,10);
yup = filter(conv(hrc, 1),1,xup); %fir1(9,0.8*1/N)


%% Agregamos ISI y ruido
% 1) Agrego ISI
ch = fir1(2, 0.4); 
ych = filter(ch,1,yup);

% 2) Del EbNo, obtengo la SNR
SNR_dB = EbNo - 10*log10(Nch) + 10*log10(k);
SNR = 10^(SNR_dB/10);
Pt = var(ych);
No = Pt/SNR;
sigma = sqrt(No/2);
ruido = sigma.*randn(length(ych),1)+1j.*sigma.*randn(length(ych),1);

rx1 = ych+ruido;

h = 1;%fir1(5, 0.9/Nch); % Anti alias por el ruido 
rx = filter(h,1,rx1);
rxch = rx(5:Nch/N:end);

%% AGC
% Mantiene la potencia de entrada del FFE constante
% para evitar que se mueva el paso optimo
target = 0.3; %VRMS
metric = std(rxch);
agc_gain = target/metric;
rx_norm = rxch(1:end).*agc_gain;

%% ECUALIZADOR FSE 
% Vamos a empezar solo con DD (slicer) sin FCR (carrier recovery)
NTAPS = 63; % Por comodidad impar
step_cma = 2^-10;%2e-3;
step_dd = 2^(-8);
i_equalizer = rx_norm(2:end); % Columna
Lrx = length(i_equalizer);

buffer = zeros(NTAPS,1); %Columna
o_filter = zeros(Lrx,1); % salida del FSE (tasa 2)
yk = zeros(fix(Lrx/N),1); % salida del FSE (tasa 1)
ak = zeros(fix(Lrx/N),1); % salida del slicer (tasa 1)
W = zeros(NTAPS,1); W((NTAPS+1)/2)=20.0;
tap_leak_gain = 1e-10;%1e-3;
force_cma_enable = 0;

error = zeros(size(ak));
coeffs = zeros( [length(ak), length(W)]);

% Arranco con CMA, luego de 1/3 de la simulacion, paso a DD
% El ultimo tercio de la simulacion la uso para procesar BER

for m=1:Lrx
    
    if (m < fix(Lrx/3))
        enable_cma=1;
        step=step_cma;
    else
        enable_cma=0;
        step=step_dd;
    end
    
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
        if enable_cma || force_cma_enable
            Ek = -yfilt*(abs(yfilt)-RCMA);
            
        else
            Ek = slicer_out-yfilt;
           
        end
        error(m/N) = Ek;
        coeffs(m/N,:) = W(:);
        
        % Emulo el sobre-muestreo del error, ejecutando el LMS
        % dentro de este if
        W = W*(1-step*tap_leak_gain) + step.*Ek.*conj(buffer);
    end
end


%% Plots de interes
%Respuesta al impulso del FFE
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
plot(f/1e9, 20*log10(abs(fftshift(fft(W/sum(W),NFFT)))));
hold all
plot(f/1e9, 20*log10(abs(fftshift(fft(h,NFFT)))));
grid on
plot(f/1e9, 20*log10(abs(fftshift(fft(hrc/sum(hrc),NFFT)))), '--k');
ylabel('Magnitude [dB]')
xlabel('Freq [GHz]')
title('Modulo de la rta. en fcia. del FFE')

figure
mse = abs(error).^2;
mse = filter(ones(1e3,1)./1e3, 1, mse);
plot(10*log10(mse))
title('MSE')
grid on
xlabel('Tiempo [samples]')
ylabel('MSE [dB]')

%% Coeffs
figure
plot(abs(coeffs))
title('Coeffs')
grid on
xlabel('Tiempo [samples]')
ylabel('Coeffs Amplitude')

%%
% Constelacion a la salida del FFE (tasa 1)
% Siempre que hagan algo adaptivo evitar las primeras muestras
A = yk(fix(2*L/3):end-100);
B = ak(fix(2*L/3):end-100);
figure
plot(real(A), imag(A),'.b'); hold all
plot(real(B), imag(B),'xr');
grid on


%% ALINEO
T0 = fix(2*L/3);
delay = finddelay(x,ak);
x2 = ak(T0+delay:end-1000);
x1 = x(T0:T0+length(x2)-1);

% SLICER
% Calcula SER simulado (no voy a calcular BER real por comodidad)
% Uds. calculen el BER como debe ser
SER= sum(x2~=x1)/length(x1);
ber_aprox = 1/k*SER
ber_teorico = berawgn(EbNo, 'qam', M)



