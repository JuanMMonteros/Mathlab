%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero
%-----------------------------------------------------------------------------%

clear; close all; clc;

%% TX Filter 

BR = 32e9;      % Baud Rate
N = 4;          % Oversampling rate
rolloff = 0.5;  % Pulse shaping rolloff
h_taps = 101;   % Pulse shaping taps

fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 consecutive samples at Tx output

% RRC
h_rrc = root_raised_cosine(BR/2, fs, rolloff, h_taps, 0);

% RRC * RRC
h_rrc_rrc = conv(h_rrc, h_rrc);

% RC
h_rc = raised_cosine(BR/2, fs, rolloff, length(h_rrc_rrc), 0);

% FFTs
NFFT = 2048;
f = (-NFFT/2:NFFT/2-1)*fs/NFFT;
H_RRC = fftshift(abs(fft(h_rrc, NFFT)));
H_RRC_RRC = H_RRC.*H_RRC;
H_RC = fftshift(abs(fft(h_rc, NFFT)));

%% Plot

% Time
n_rrc_v = -(h_taps-1)/2:(h_taps-1)/2;
n_rc_v = -(length(h_rrc_rrc)-1)/2:(length(h_rrc_rrc)-1)/2;
figure; 
plot(n_rrc_v, h_rrc, '--k', 'Linewidth', 1);
hold on; grid on;
plot(n_rc_v,h_rrc_rrc, 'r', 'Linewidth', 2);
plot(n_rc_v,h_rc, '--b', 'Linewidth', 1.5);
title('hrrc*hrrc = hrc'); 
xlabel('Samples'); ylabel('Amplitude')
legend('hrrc','hrrc*hrrc','hrc')
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

% Freq
figure; 
plot(f/1e9, H_RRC, '--k', 'Linewidth', 1);
hold on; grid on;
plot(f/1e9, H_RRC_RRC, 'r', 'Linewidth', 2);
plot(f/1e9, H_RC, '--b', 'Linewidth', 1.5);
title('Hrrc.Hrrc = Hrc'); 
xlabel('Freq [GHz]'); ylabel('Amplitude')
legend('Hrrc','Hrrc.Hrrc','Hrc')
set(gcf, 'Position', [50 50 500 500],'Color', 'w');
