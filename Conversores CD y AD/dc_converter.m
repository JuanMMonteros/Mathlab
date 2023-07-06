%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% This script simulates a D/C convertor

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

%% Parameters

fs = 1e6;               % Sampling frequency [Hz]
Ts = 1/fs;              % Sampling time [s]

OVR = 8;                % Oversampling factor (fs_ch/fs)
fs_ch = fs*OVR;         % Channel sampling frequency [Hz]
 
n_samples = 1e2;        % Sim length
f0 = 1e5;               % Tone freq

fil_sel = 0;            % Filter selector. 0:ZOH; 1:FIR1
h_taps = 151;           % Filter taps 
h_fc = 0.4*fs;          % Filter cut frequency

fz = 15;                % Plots Font Size

%% Processing

% Time vector
t_v = (0:n_samples-1)*Ts;
t_up_v = (0:n_samples*OVR-1)*1/fs_ch;

% Signal in digital domain
x_d_v = 2*cos(2*pi*f0*t_v); 

% Ideal signal in continuous domain
x_ci_v = 2*cos(2*pi*f0*t_up_v);

% Upsampling
x_up_v = upsample(x_d_v, OVR);

% Reconstruction filter
switch fil_sel

    case 0
        h_v = ones(1,OVR);
        h_delay = 0; % Filter delay
    case 1
        h_v = fir1(h_taps-1, h_fc/(fs_ch/2));
        h_delay =(h_taps-1)/2; % Filter delay
    otherwise
        error('Wrong configuration')
end
  
h_v = OVR * h_v / sum(h_v); % Gain to keep constant the signal amplitude

% Filtering
x_c_v = filter(h_v, 1, [x_up_v, zeros(1,h_delay)]);
x_c_v = x_c_v(h_delay+1:end);

% FFT variables
NFFT = 16 * n_samples;
f_v = (-NFFT/2:NFFT/2-1)*fs/NFFT;
f_up_v = (-NFFT*OVR/2:NFFT*OVR/2-1)*fs_ch/(NFFT*OVR);

% FFTs
X_d_v = fftshift(abs(fft(x_d_v, NFFT))/n_samples);
X_up_v = fftshift(abs(fft(x_up_v, NFFT*OVR))/(n_samples*OVR));
X_c_v = fftshift(abs(fft(x_c_v, NFFT*OVR))/(n_samples*OVR));
X_ci_v = fftshift(abs(fft(x_ci_v, NFFT*OVR))/(n_samples*OVR));
H_v = fftshift(abs(fft(h_v, NFFT*OVR)));
H_v = H_v/max(H_v);

%% Plots 
% Time
figure;
p = plot(t_v/1e-6, x_d_v ,'ob','Linewidth',1);
p.MarkerFaceColor = p.Color;
p.MarkerEdgeColor = p.Color;
hold on;
plot(t_up_v/1e-6, x_up_v ,'xk','Linewidth',1);
plot(t_up_v/1e-6, x_c_v ,'-r','Linewidth',1.5);
plot(t_up_v/1e-6, x_ci_v ,'--c','Linewidth',1);
tit = sprintf('D/C converter');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Time [us]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'D','UP','C','C-ideal'}, 'Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

% Freqfigure;
figure
plot(f_v/1e6, X_d_v ,'-b','Linewidth',2);
hold on;
plot(f_up_v/1e6, X_up_v ,'-.k','Linewidth',2);
plot(f_up_v/1e6, X_c_v ,'--r','Linewidth',2);
plot(f_up_v/1e6, X_ci_v ,'-.c','Linewidth',1);
plot(f_up_v/1e6, H_v ,'-.m','Linewidth',1);
tit = sprintf('D/C converter');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Frequency [MHz]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'D','UP','C','C-ideal','H'}, 'Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [550 50 800 500],'Color', 'w');
