%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% This script shows the FFT resolution of a tone

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

%% Parameters

en_plots = 1;       % 1: ON| 0: OFF
en_second_tone = 0; % 1: ON| 0: OFF

fs = 2e9;           % Sampling frequency [Hz]
ts = 1/fs;          % Sampling time [s]

f0 = 100e6;         % Tone frequency [Hz]
T = 1/f0;           % Tone period [s]
Ps_dbm = 10;        % Signal power [dBm] 
t_meas = 100e-9;    % Signal length [s]

fft_zp = 8;         % FFT zero pading factor

fz = 15;            % Plots Font Size

%% Processing

n_samples = round(t_meas*fs);

% Tone resolution
delta_f = 1/t_meas;

% Time vector
t_v = (0:n_samples-1)*ts;

% dBm to W
Ps = 1e-3*10^(Ps_dbm/10);  % Signal power [W]
   
% Signal
x_v = sqrt(Ps)*exp(1j*2*pi*f0*t_v);

if en_second_tone
    x_v = x_v + sqrt(Ps)*exp(1j*2*pi*(f0+2*delta_f)*t_v);
end

% FFT
NFFT = n_samples * fft_zp;
f_v = (0:NFFT-1)*fs/NFFT;
X_v = fft(x_v, NFFT)/n_samples;


%% Print
fprintf('\n --------------------------------\n');
fprintf(' - Resolution theory: %.2f [MHz]\n', delta_f);
fprintf(' --------------------------------\n');
   
%% Plots

% Signal PSD

if en_plots

    figure;
    subplot(2,1,1)
    %plot(f_v/1e6,10*log10(abs(X_v)/1e-3),'-r','Linewidth',2); hold on;
    plot(f_v/1e6,(abs(X_v)/1e-3),'-r','Linewidth',2); hold on;
    ylims = ylim;
    plot([f0-delta_f,f0+delta_f]/1e6,[ylims(2)/2,ylims(2)/2],'--m','Linewidth',2)
    tit = sprintf('Tone plot. \n T = %.2f[ns]. f0 = %.2f[MHz] \n Tmeas = %.2f[us]. Res = %.2f[MHz]',...
                                                                T/1e-9, f0/1e6,t_meas/1e-6, delta_f/1e6);
    title(tit,'Interpreter','latex','FontSize', fz);
    legend({'FFT','2/Tmeas'}, 'Interpreter','latex','FontSize', fz-2);
    xlabel('Frequency [MHz]', 'Interpreter','latex','FontSize', fz);
    ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
    grid on; xlim([f0/1e6/2,f0/1e6*2])
    subplot(2,1,2)
    plot(t_v/1e-9,real(x_v),'-r','Linewidth',2); hold on;
    plot(t_v/1e-9,imag(x_v),'-b','Linewidth',2);
    legend({'R','I'}, 'Interpreter','latex','FontSize', fz-2);
    xlabel('Time [ns]', 'Interpreter','latex','FontSize', fz);
    ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
    grid on; 
    set(gcf, 'Position', [50 50 700 900],'Color', 'w');
    
end