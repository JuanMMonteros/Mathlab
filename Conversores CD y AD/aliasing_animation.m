%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% This script shows aliasing animation

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

%% Parameters

fs = 1e6;               % Sampling frequency [Hz]
Ts = 1/fs;              % Sampling time [s]

OVR = 5;                % Oversampling factor (fs_ch/fs)
fs_ch = fs*OVR;         % Channel sampling frequency [Hz]

Ps_dbm = 20;            % Signal power [dBm] 
n_samples = 1e4;        % Sim length
max_freq = 0.45*fs_ch;  % Max tone freq

n_exp = 30;             % Number of experiments

fz = 15;                % Plots Font Size

%% Processing

% Tone freq vector to sweep
f0_v = (0:n_exp-1)*max_freq/(n_exp-1);   

% Time vector
t_up_v = (0:n_samples*OVR-1)*1/fs_ch;
t_v = (0:n_samples-1)*1/Ts;

% dBm to W
Ps = 1e-3*10^(Ps_dbm/10);  % Signal power [W]

% FFT variables
NFFT = 16 * n_samples;
f_v = (0:NFFT-1)*fs/NFFT;
f_up_v = (0:NFFT*OVR-1)*fs_ch/(NFFT*OVR);

for idx = 1:n_exp
    
    % Print every 10 experiments 
    if mod(idx,10) == 0
       fprintf('- Progress:  %.0f %% \n', 100*idx/n_exp);
    end
    
    % Tone frequency selector
    f0 = f0_v(idx);
   
    % Signals
    x1_up_v = sqrt(2*Ps)*cos(2*pi*f0*t_up_v);   % A = sqrt(2*P)
    x1_v = x1_up_v(1:OVR:end);
    
    x2_up_v = sqrt(Ps)*exp(1j*2*pi*f0*t_up_v);  % A = sqrt(P)
    x2_v = x2_up_v(1:OVR:end);

    % FFTs
    X1_v = (abs(fft(x1_v, NFFT))/n_samples).^2;
    X1_up_v = (abs(fft(x1_up_v, NFFT*OVR))/(n_samples*OVR)).^2;
    
    X2_v = (abs(fft(x2_v, NFFT))/n_samples).^2;
    X2_up_v = (abs(fft(x2_up_v, NFFT*OVR))/(n_samples*OVR)).^2;
    
    % -- FFTs advance plot --
    % X1
    figure(1); clf;
    plot((f_up_v-fs_ch/2)/1e6,fftshift(X1_up_v),'-b','Linewidth',2)
    hold on;
    plot((f_v-fs/2)/1e6,fftshift(X1_v),'--r','Linewidth',2)
    tit = sprintf('FFT advance. f0 =  %.1e, fs =  %.1e', f0, fs);
    title(tit,'Interpreter','latex','FontSize', fz);
    xlabel('Frequency [MHz]', 'Interpreter','latex','FontSize', fz);
    ylabel('Power [W]', 'Interpreter','latex','FontSize', fz);
    grid on
    ylim([0,1.6*Ps])
    set(gcf, 'Position', [50 50 500 500],'Color', 'w');
    
    % X2
    figure(2); clf;
    plot((f_up_v-fs_ch/2)/1e6,fftshift(X2_up_v),'-b','Linewidth',2)
    hold on;
    plot((f_v-fs/2)/1e6,fftshift(X2_v),'--r','Linewidth',2)
    tit = sprintf('FFT advance. f0 =  %.1e, fs =  %.1e', f0, fs);
    title(tit,'Interpreter','latex','FontSize', fz);
    xlabel('Frequency [MHz]', 'Interpreter','latex','FontSize', fz);
    ylabel('Power [W]', 'Interpreter','latex','FontSize', fz);
    grid on
    ylim([0,1.6*Ps])
    set(gcf, 'Position', [600 50 500 500],'Color', 'w');
    
    pause(.1) % Pause .1 [sec]
    
end