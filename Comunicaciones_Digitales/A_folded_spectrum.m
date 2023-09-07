clear all
close all

%Basic TX QPSK 
BR = 32e9; %Bd
L = 10000; % Simulati on Length
N = 2; % Oversampling rate
fs = N*BR; % Sampling rate to emulate analog domain
T = 1/BR; % Time interval between two consecutive symbols
Ts = 1/fs; % Time between 2 conseutive samples at Tx output
EbNo = 6; %dB
k = 2; %k=log2(M)

% Two symbols generation (+1,-1)
xi = 2*randi([0,1],L,1)-1;
xq = 2*randi([0,1],L,1)-1;
x=xi+j.*xq;

% Upsampling to change sampling rate
xup = upsample(x,N);

% Filter to interpolate the signal
g = rcosine(BR,fs,'sqrt', 0.5,10);%coseno_realzado(BR, fs, 0.8, 100);
yup = filter(g,1,xup);

%% Agrego ISI y ruido
%b = fir1(2,1.2/N); % Si el orden es impar hay que re-ajustar la fase de muestreo
b = 1;
ych = filter(b,1,yup);

SNR_dB = EbNo - 10*log10(N) + 10*log10(k);
SNR = 10^(SNR_dB/10);
Pt = var(yup);
No = Pt/SNR;
sigma = sqrt(No/2);
ruido = sigma.*randn(length(yup),1)+1j.*sigma.*randn(length(yup),1);

rx = ych+ruido;

%% MATCHED FILTER h*(-t)
matched_filter = conj(conv(b,g)); matched_filter=matched_filter(end:-1:1);
ymf = filter(matched_filter,1,rx);
rx_down = ymf(1:N:end);

figure; plot(real(rx_down),'.');

scatterplot(rx_down(100:end));
title('Constelacion a la salida del MF');

% Respuesta equivalente del SMF h(t) x h*(-t)|t=kT  = rho(k)
% El MF acentua el ISI!
figure
h = conv(b,g);
ht = conv(h,matched_filter);
rho_k = downsample(ht,N);
stem(rho_k);
title('Respuesta equivalente del SMF')
xlabel('Samples'); ylabel('Rho(k) Amplitude');

%%
NFFT = 2048;
NOS = 1;

fft_channel = fft(upsample(h,NOS),NFFT);
fft_channel_mf = fft(upsample(ht,NOS),NFFT);
fft_folded_spectrum = fft( upsample(rho_k,NOS*N), NFFT);
f = -fs*NOS/2:NOS*fs/NFFT:fs*NOS/2-NOS*fs/NFFT;

figure
freqr = abs(fftshift(fft_channel)); freqr=freqr/freqr(NFFT/2);
plot(f/1e9,freqr)
hold all
freqr = abs(fftshift(fft_channel_mf)); freqr=freqr/freqr(NFFT/2);
plot(f/1e9,freqr)
freqr = abs(fftshift(fft_folded_spectrum)); freqr=freqr/freqr(NFFT/2);
plot(f/1e9,freqr, 'LineWidth', 3)
grid on
legend('FFT {h(t)}', 'FFT {h(t) x h*(-t)}', 'Folded Spectrum')
title('Frequency Response of Folded Spectrum')
ylabel('Amplitude [times]')
xlabel('Frequency [GHz]')

figure
pwelch(rx_down, hanning(1024), 0, 4096)

