%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% This script simulates a C/D convertor

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

%% Parameters

fs = 5e9;               % Sampling frequency [Hz]
Ts = 1/fs;              % Sampling time [s]

OVR = 4;                % Oversampling factor (fs_ch/fs)
fs_ch = fs*OVR;         % Channel sampling frequency [Hz]
 
n_samples = 1e2;        % Sim length in digital domain
omega_0 = 2*pi*100e6;   % Tone frequency

h_taps = 151;           % Filter taps 
h_fc = 0.4*fs;          % Filter cut frequency

fz = 15;                % Plots Font Size

%% Processing

% Time vector
t_up_v = (0:n_samples*OVR-1)*1/fs_ch;
t_v = (0:n_samples-1)*Ts;

% Ideal signal in digital domain
tone_dig = 2*cos(2*pi*omega_0*t_v); 

% Signal in contiuous domain
tone_cont = 2*cos(2*pi*omega_0*t_up_v);  

% AntiAlias filter (AAF)
filter_aa = fir1(h_taps-1, h_fc/(fs_ch/2));
h_delay =(h_taps-1)/2; % Filter delay

% Filtering
tone_cont_f = filter(filter_aa, 1, [tone_cont, zeros(1,h_delay)]);
tone_cont_f = tone_cont_f(h_delay+1:end);

% -- Downsampling --
% Option 1
%x_d_v = downsample(x_cf_v, OVR);  
% Option 2
tone_cf_dig = tone_cont_f(1:OVR:end);   

% FFT variables
NFFT = 16 * n_samples;
f_v = (-NFFT/2:NFFT/2-1)*fs/NFFT;
f_up_v = (-NFFT*OVR/2:NFFT*OVR/2-1)*fs_ch/(NFFT*OVR);

% FFTs
X_d_v = fftshift(abs(fft(tone_cf_dig, NFFT))/n_samples);
X_di_v = fftshift(abs(fft(tone_dig, NFFT))/(n_samples));
X_cf_v = fftshift(abs(fft(tone_cont_f, NFFT*OVR))/(n_samples*OVR));
X_c_v = fftshift(abs(fft(tone_cont_f, NFFT*OVR))/(n_samples*OVR));
H_v = fftshift(abs(fft(filter_aa, NFFT*OVR)));

%% Plots 

% Time
figure;
p = plot(t_up_v/1e-6, tone_cont ,'-r','Linewidth',1);
p.MarkerFaceColor = p.Color;
p.MarkerEdgeColor = p.Color;
hold on;
plot(t_up_v/1e-6, tone_cont_f ,'--k','Linewidth',1);
plot(t_v/1e-6, tone_cf_dig ,'ob','Linewidth',1.5);
plot(t_v/1e-6, tone_dig ,'--c','Linewidth',1);
tit = sprintf('C/D converter');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Time [us]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'C','C-fil','D','D-ideal'}, 'Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

% Freq
figure
plot(f_up_v/1e6, X_c_v ,'-r','Linewidth',2);
hold on;
plot(f_up_v/1e6, X_cf_v ,'--k','Linewidth',2);
plot(f_v/1e6, X_d_v ,'--b','Linewidth',2);
plot(f_v/1e6, X_di_v ,'-.c','Linewidth',1);
plot(f_up_v/1e6, H_v ,'-.m','Linewidth',1);
tit = sprintf('C/D converter');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Frequency [MHz]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'C','Fil','D','D-ideal','H'}, 'Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [550 50 800 500],'Color', 'w');

