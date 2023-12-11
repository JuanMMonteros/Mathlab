%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero
%-----------------------------------------------------------------------------%

close all; clear

%% Basic TX QPSK 
L = 10000;      % Simulation Length
BR = 32e9;      % Baud Rate
N = 4;          % Oversampling rate
rolloff = 0.1;  % Pulse shaping rolloff
h_taps = 101;   % Pulse shaping taps

fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 consecutive samples at Tx output

% Two symbols generation (+1,-1) for QPSK
xi = 2*randi([0,1],L,1)-1; % Real
xq = 2*randi([0,1],L,1)-1; % Imag
x = xi + 1j*xq;

% Upsampling to change sampling rate
xup = N * upsample(x,N);

% Filter to interpolate the signal
h = raised_cosine(BR/2, fs, rolloff, h_taps, 0);
h_delay = (h_taps-1)/2;
yup = filter(h,1,[xup; zeros(h_taps, 1)]);
yup = yup(1+h_delay:end);

%%
% Plot of Tx output and original symbols
t_symbols = [0:L-1].*T; %Timeline for symbols before upsampling 
t_samples = [0:length(yup)-1].*Ts;
LPLOT = 30*N; % Total symbols to plot
figure

subplot(2,1,1)
plot(t_samples, real(yup),'-x', 'Linewidth',2)
hold all
stem(t_symbols, real(x), 'Linewidth',2)
ylabel('Amplitude')
xlabel('Time in seconds')
title('Tx Real part')
grid on
xlim([0, t_samples(LPLOT)])

subplot(2,1,2)
plot(t_samples, imag(yup),'-x', 'Linewidth',2)
hold all
stem(t_symbols, imag(x), 'Linewidth',2)
ylabel('Amplitude')
xlabel('Time in seconds')
title('Tx Imag part')
grid on
xlim([0, t_samples(LPLOT)])
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

%% Plot filter design
NFFT = 1024*8;

figure; 

subplot(2,1,1)
stem(h, 'r', 'Linewidth',2)
title('h')
xlabel('Samples')
ylabel('Amplitude')
xlim([0,length(h)]); grid on;

subplot(2,1,2)
H_abs = abs(fft(h,NFFT));
f = (-NFFT/2:NFFT/2-1)*fs/NFFT;
plot(f/1E9, fftshift(H_abs), '-r', 'Linewidth',2);
hold on;
plot([BR/2/1E9 BR/2/1E9], ylim, '--k', 'Linewidth',1)
plot([BR/2/1E9 BR/2/1E9]*(1+rolloff), ylim, '--m', 'Linewidth',1)
title('|H|')
xlabel('Freq [GHz]')
ylabel('Amplitude')
grid on;

legend({'|H|','BR/2','BR/2*(1+Beta)'},'Location','nw')

%% PSD AT TX OUTPUT
NFFT = 1024*8;
WELCH_OVERLAP = 0.*NFFT;

figure
[Pxx, f] = pwelch(x, hanning(NFFT/2), WELCH_OVERLAP, NFFT, BR, 'centered');
Px_dB= 10*log10(Pxx);
Px_dB = Px_dB - max(Px_dB); %Normalize Px_dB to get 0dB at f=0, only for plot purposes
plot(f/1e9, (Px_dB),'-b', 'Linewidth',1)
hold on

[Pxx, f] = pwelch(x, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs, 'centered');
Px_dB= 10*log10(Pxx);
Px_dB = Px_dB - max(Px_dB); %Normalize Px_dB to get 0dB at f=0, only for plot purposes
plot(f/1e9, (Px_dB),'-.m', 'Linewidth',1)
hold on

[Pxx, f] = pwelch(yup, hanning(NFFT/4), WELCH_OVERLAP, NFFT, fs, 'centered');
Px_dB= 10*log10(Pxx);
Px_dB = Px_dB - max(Px_dB); %Normalize Px_dB to get 0dB at f=0, only for plot purposes
plot(f/1e9, Px_dB,'-r', 'Linewidth',1)
grid on
xlabel('Frequency [GHz]')
ylabel('PSD Magnitude [dB(V^2)/Hz]')

hold on
H_abs = abs(fft(h,NFFT));
H_dB = 20.*log10(H_abs/H_abs(1)); % It is a filter, is 20log10, PSD is power ant it is 10log10
plot(f/1E9, fftshift(H_dB), '--k', 'Linewidth',2)

legend({'PSD symbols','PSD symbols up','PSD symbols fil','20*log10(H)'},'Location','s')
title('PSDs')

set(gcf, 'Position', [50 50 1000 500],'Color', 'w');

%% Diagrama de ojo
eyediagram(yup(500:5000),N*2)
set(gcf, 'Position', [50 50 500 500],'Color', 'w');