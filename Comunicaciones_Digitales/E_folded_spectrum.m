%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero
%-----------------------------------------------------------------------------%

clear; close all; clc;

 

BR = 32e9;      % Baud Rate
N = 4;          % Oversampling rate
rolloff = 0.5;  % Pulse shaping rolloff
h_taps = 101;   % Pulse shaping taps

fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 consecutive samples at Tx output

%% TX Filter

h_tx = root_raised_cosine(BR/2, fs, rolloff, h_taps, 0);

%% CH filter

h_ch = [0 0 0 0 1 0 0 0 0]; % Impulse
%h_ch = fir1(1, 0.9);

h_tx_ch = conv(h_tx, h_ch);

%% RX filter 

h_rx = conj(h_tx_ch(end:-1:1));

%% Total filter

h_tot = conv(h_tx_ch, h_rx);

h_tot_dw = N*h_tot(1:N:end);
% Si h_tot_dw es un impulso h_tot es un pulso de Nyquist 

%% Spectrums
NFFT = 2048;

f = (-NFFT/2:NFFT/2-1)*fs/NFFT;

H_tx = fftshift(abs(fft(h_tx, NFFT)));
H_ch = fftshift(abs(fft(h_ch, NFFT)));
H_tx_ch = fftshift(abs(fft(h_tx_ch, NFFT)));
H_rx = fftshift(abs(fft(h_rx, NFFT)));
H_tot = fftshift(abs(fft(h_tot, NFFT)));

f_dw = (-NFFT/2:NFFT/2-1)*BR/NFFT;
H_tot_dw = fftshift(abs(fft(h_tot_dw, NFFT)));
% Si H_tot_dw es plana entonces h_tot es un pulso de Nyquist 

%% Plots

% Time
figure; 

subplot(2,3,1); 
stem(h_tx); xlim([0,length(h_tx)]); grid on
title('Tx'); xlabel('Samples'); ylabel('Amplitude')

subplot(2,3,2);
stem(h_ch); xlim([0,length(h_ch)]); grid on
title('Ch'); xlabel('Samples'); ylabel('Amplitude')

subplot(2,3,3);
stem(h_tx_ch); xlim([0,length(h_tx_ch)]);grid on
title('Tx+Ch'); xlabel('Samples'); ylabel('Amplitude')

subplot(2,3,4);
stem(h_rx); xlim([0,length(h_rx)]); grid on
title('Rx'); xlabel('Samples'); ylabel('Amplitude')

subplot(2,3,5);
stem(h_tot); xlim([0,length(h_tot)]); grid on
title('Tx+Ch+Rx'); xlabel('Samples'); ylabel('Amplitude')

subplot(2,3,6);
stem(h_tot_dw); xlim([0,length(h_tot_dw)]); grid on
title('(Tx+Ch+Rx) DW'); xlabel('Samples'); ylabel('Amplitude')

set(gcf, 'Position', [50 50 1000 600],'Color', 'w');

% Freq
figure; 

plot(f/1e9,abs(H_tot),'-r', 'Linewidth',1);         % Espectro original
hold on
plot((f+1*BR)/1e9,abs(H_tot),'-r', 'Linewidth',1)   % Replica en BR
plot((f+2*BR)/1e9,abs(H_tot),'-r', 'Linewidth',1)   % Replica en 2BR
plot((f+3*BR)/1e9,abs(H_tot),'-r', 'Linewidth',1)   % Replica en 3BR
plot((f-1*BR)/1e9,abs(H_tot),'-r', 'Linewidth',1)  % Replica en -BR
plot((f-2*BR)/1e9,abs(H_tot),'-r', 'Linewidth',1)  % Replica en -2BR
plot((f-3*BR)/1e9,abs(H_tot),'-r', 'Linewidth',1)  % Replica en -3BR
plot(f_dw/1e9, H_tot_dw, '-k', 'Linewidth', 2);    % Espectro plegado
plot((f_dw+1*BR)/1e9, H_tot_dw, '-k', 'Linewidth', 2);    % Espectro plegado
plot((f_dw+2*BR)/1e9, H_tot_dw, '-k', 'Linewidth', 2);    % Espectro plegado
plot((f_dw+3*BR)/1e9, H_tot_dw, '-k', 'Linewidth', 2);    % Espectro plegado
plot((f_dw-1*BR)/1e9, H_tot_dw, '-k', 'Linewidth', 2);    % Espectro plegado
plot((f_dw-2*BR)/1e9, H_tot_dw, '-k', 'Linewidth', 2);    % Espectro plegado
plot((f_dw-3*BR)/1e9, H_tot_dw, '-k', 'Linewidth', 2);    % Espectro plegado
title('Folded spectrum'); grid on;
xlabel('Freq [GHz]'); ylabel('Amplitude')
set(gcf, 'Position', [50 50 500 500],'Color', 'w');
set(gcf, 'Position', [50 50 600 600],'Color', 'w');

