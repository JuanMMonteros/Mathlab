%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero
%-----------------------------------------------------------------------------%

clear 
close all
clc;

%% Basic TX QPSK 
L = 10000;      % Simulation Length
BR = 32e9;      % Baud Rate
N = 4;          % Oversampling rate
rolloff = 0.5;  % Pulse shaping rolloff
h_taps = 101;   % Pulse shaping taps

fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 consecutive samples at Tx output

% Two symbols generation (+1,-1) for QPSK
xi = 2*randi([0,1],L,1)-1; % Real
xq = 2*randi([0,1],L,1)-1; % Imag
x = xi + 1j*xq;

% Upsampling to change sampling rate
xup = upsample(x,N);

% Filter to interpolate the signal
h = raised_cosine(BR/2, fs, rolloff, h_taps, 0);
h_delay = (h_taps-1)/2;
yup = filter(h,1,[xup; zeros(h_taps, 1)]);
yup = yup(1+h_delay:end);

%% Let's transform yup to analytic signal
% yup is a real signal with BW=(1+rolloff)*BR/2
% yup must be rotated in frequency domain until the whole spectrum is
% placed at positive frequencies
t = (0:length(yup)-1).'.*Ts; %Note the operator .' to transpose
f0 = (1+rolloff)*BR/2* (1.2); %20% extra shift just in case
yup_analytic = yup.*exp(1j*2*pi*f0.*t);

%% Plots

figure
NFFT = 1024;
WELCH_OVERLAP = 0;

[Pxx, f] = pwelch(yup, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs);
Px_dB= 10*log10(Pxx);
Px_dB = Px_dB - Px_dB(1); %Normalize Px_dB to get 0dB at f=0, only for plot purposes
plot(f, Px_dB,'-b', 'Linewidth',1)
grid on
xlabel('Frequency [Hz]')
ylabel('PSD Magnitude [dB(V^2)/Hz]')
title('PSD Original signal')
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

figure
[Pxx, f] = pwelch(yup_analytic, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs);
Px_dB= 10*log10(Pxx);
Px_dB = Px_dB - Px_dB(1); %Normalize Px_dB to get 0dB at f=0, only for plot purposes
plot(f, Px_dB,'-b', 'Linewidth',1)
grid on
xlabel('Frequency [Hz]')
ylabel('PSD Magnitude [dB(V^2)/Hz]')
title('PSD Analytic signal')
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

%% Split analytic signal in real and imaginary part
y_real = real(yup_analytic);
y_imag = imag(yup_analytic);

%% Model hilber transform
analytic_signal = hilbert(y_real);
figure
plot(imag(analytic_signal), '-b');
hold all
plot(y_imag, 'xr--');
title('Hilbert transform')
xlabel('Samples')
ylabel('Amplitude')
legend('Hilbert output','Original')
xlim([0,100]); grid on;
set(gcf, 'Position', [50 50 500 500],'Color', 'w');
