%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% The nomenclature "_v" at the end of a variable indicates that is a vector
% The nomenclature "_m" at the end of a variable indicates that is a matrix
% The nomenclature "_s" at the end of a variable indicates that is a struct
% The nomenclature "_c" at the end of a variable indicates that is a cell

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

%%

% WARNING: be careful with the PC volume when running this script

% This script reads an audio signal with an interference in 1KHz and try to
% eliminate it by filtering the signal

% Read signal
[s_v,fs] = audioread('signal_with_interference.wav');
n_samples = length(s_v);
time_len = n_samples*1/fs;

% FFT
NFFT = 8*n_samples;
S_v = fft(s_v, NFFT)/n_samples;
S_db_v = 10*log10(abs(S_v).^2);

% Time and frequency vectors
t_v = (0:n_samples-1)*1/fs;
f_v = (0:NFFT-1)*fs/NFFT;

% Plot input signals in time
fz = 15; % Fontsize
figure;
plot(t_v,s_v,'-r','Linewidth',2)
title('Input in time','Interpreter','latex','FontSize', fz);
xlabel('Time [s]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
grid on
set(gcf, 'Position', [50 50 400 400],'Color', 'w');

% Signals in freq
figure;
plot((f_v-fs/2)/1e3,fftshift(S_db_v),'-r','Linewidth',2)
title('Input in frequency','Interpreter','latex','FontSize', fz);
xlabel('Frequency [KHz]', 'Interpreter','latex','FontSize', fz);
ylabel('Module [dB]', 'Interpreter','latex','FontSize', fz);
grid on; 
set(gcf, 'Position', [500 50 600 600],'Color', 'w');

% Play the input signal and pause the simulation until it ends
sound(s_v, fs);
pause(time_len);

%% Filter design

% Filter desing: make a Notch filter in 1KHz
f_notch = 1e3;
f0 = f_notch/(fs/2);
[b_notch, a_notch] = iirnotch(f_notch/(fs/2), f0/50);
H_v = freqz(b_notch,a_notch,NFFT,'whole',fs);
H_db_v = 20*log10(abs(H_v));

% Filter response
figure;
subplot(2,1,1)
plot((f_v-fs/2)/1e3,fftshift(H_db_v),'-r','Linewidth',1.5)
title('Filter response','Interpreter','latex','FontSize', fz);
xlabel('Frequency [KHz]', 'Interpreter','latex','FontSize', fz);
ylabel('Module [dB]', 'Interpreter','latex','FontSize', fz);
grid on; 
subplot(2,1,2)
plot((f_v-fs/2)/1e3,fftshift(angle(H_v)),'-r','Linewidth',1.5)
xlabel('Frequency [KHz]', 'Interpreter','latex','FontSize', fz);
ylabel('Angle', 'Interpreter','latex','FontSize', fz);
grid on;
set(gcf, 'Position', [500 50 600 600],'Color', 'w');

%% Filtering

s_fil_v = filter(b_notch, a_notch, s_v); 

% FFT
S_fil_v = fft(s_fil_v, NFFT)/n_samples;
S_fil_db_v = 10*log10(abs(S_fil_v).^2);

% Plot input signals in time
fz = 15; % Fontsize
figure;
plot(t_v,s_fil_v,'-b','Linewidth',2)
title('Output in time','Interpreter','latex','FontSize', fz);
xlabel('Time [s]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
grid on
set(gcf, 'Position', [50 50 400 400],'Color', 'w');

% Signals in freq
figure;
plot((f_v-fs/2)/1e3,fftshift(S_db_v),'-r','Linewidth',2)
hold on;
plot((f_v-fs/2)/1e3,fftshift(S_fil_db_v),'-b','Linewidth',2)
title('Signals in frequency','Interpreter','latex','FontSize', fz);
legend({'Input','Output'}, 'Interpreter','latex','FontSize', fz-2);
xlabel('Frequency [KHz]', 'Interpreter','latex','FontSize', fz);
ylabel('Module [dB]', 'Interpreter','latex','FontSize', fz);
grid on; 
set(gcf, 'Position', [500 50 600 600],'Color', 'w');

% Play the input signal and pause the simulation until it ends
sound(s_fil_v, fs);
pause(time_len);