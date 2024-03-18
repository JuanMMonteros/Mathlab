%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% This script shows how to measure the SNR of a tone using a PSD

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

%% Parameters

en_plots = 1;       % 1: ON| 0: OFF

fs = 2e9;           % Sampling frequency [Hz]
ts = 1/fs;          % Sampling time [s]

f0 = 100e6;         % Tone frequency [Hz]
Ps_dbm = 10;        % Signal power [dBm] 
t_meas = 1e-6;      % Signal length [s]
snr_db = 15;        % Signal SNR [dB]

fft_zp = 8;         % FFT zero pading factor

n_exp = 500;        % Number of experiments

fz = 15;            % Plots Font Size

%% Processing

n_samples = round(t_meas*fs);

% Time vector
t_v = (0:n_samples-1)*ts;

% dBm to W
Ps = 1e-3*10^(Ps_dbm/10);  % Signal power [W]

snr = 10^(snr_db/10);

% Noise power: SNR = Es/En = (Ps*t_meas)/(pn/BW);
% Es: Signal energy = Ps * t_meas
% En: Noise energy = Pn/BW
% BW = signal bandwidth. For complex signals BW=fs
% We can compute the signal power as: Pn = SNR/(Ps*t_meas)*BW

Pn = (Ps*t_meas)/snr*fs;

% FFT
NFFT = n_samples * fft_zp;
f_v = (0:NFFT-1)*fs/NFFT;

% PSD
psd_v = 0;

for idx = 1:n_exp
    
    % Print every 10 experiments 
    if mod(idx,10) == 0
       fprintf('- Progress:  %.0f %% \n', 100*idx/n_exp);
    end
   
    % Signal
    x_ideal_v = sqrt(Ps)*exp(1j*2*pi*f0*t_v);

    % Noise
    n_v = sqrt(Pn/2)*(randn(1,n_samples)+1j*randn(1,n_samples));

    % Signal + noise
    x_v = x_ideal_v + n_v;

    % FFT
    X_v = fft(x_v, NFFT)/n_samples;

    % PSD
    psd_v = psd_v + abs(X_v).^2/n_exp;

    % PSD advance plot every 10 exp
    if en_plots && mod(idx,10)==0
        figure(1); clf;
        plot(f_v/1e6,psd_v,'-b','Linewidth',2)
        tit = sprintf('PSD advance. Exp %d of %d', idx, n_exp);
        title(tit,'Interpreter','latex','FontSize', fz);
        xlabel('Frequency [MHz]', 'Interpreter','latex','FontSize', fz);
        ylabel('Power [W]', 'Interpreter','latex','FontSize', fz);
        grid on
        ylim([0,1.2*Ps])
        set(gcf, 'Position', [50 50 700 500],'Color', 'w');
        pause(.1)
    end
    
end

% SNR computation
[psd_max, max_idx] = max(psd_v);
if max_idx < n_samples/2
    Pn_est = mean(psd_v(round(NFFT/2):end));
else
    Pn_est = mean(psd_v(1:round(NFFT/2)));
end
Ps_est = psd_max - Pn_est;
snr_est = Ps_est/Pn_est;
snr_est_db = 10*log10(snr_est);

%% Print
fprintf('\n --------------------------------\n');
fprintf(' - Ps theory: %.2f [dBm]\n',Ps_dbm);
fprintf(' - Ps estimated: %.2f [dBm]\n',10*log10(Ps_est/1e-3));
fprintf(' - Pn theory: %.2f [dBm]\n',10*log10(Pn/n_samples/1e-3));
fprintf(' - Pn estimated: %.2f [dBm]\n',10*log10(Pn_est/1e-3));
fprintf(' - SNR theory: %.2f [dB]\n',snr_db);
fprintf(' - SNR estimated: %.2f [dB]\n',snr_est_db);
fprintf(' - Estimation error: %.2f [dB]\n',snr_db-snr_est_db);
fprintf(' --------------------------------\n');
   
%% Plots

% Signal PSD

if en_plots

    figure;
    subplot(2,1,1)
    plot(f_v/1e6,10*log10(abs(X_v)/1e-3),'-r','Linewidth',2); hold on;
    title('1 experiment plot','Interpreter','latex','FontSize', fz);
    xlabel('Frequency [MHz]', 'Interpreter','latex','FontSize', fz);
    ylabel('Amplitude [dBm]', 'Interpreter','latex','FontSize', fz);
    grid on; 
    subplot(2,1,2)
    plot(t_v/1e-6,real(x_v),'-r','Linewidth',2); hold on;
    plot(t_v/1e-6,imag(x_v),'-b','Linewidth',2);
    legend({'R','I'}, 'Interpreter','latex','FontSize', fz-2);
    xlabel('Time [us]', 'Interpreter','latex','FontSize', fz);
    ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
    grid on; 
    set(gcf, 'Position', [50 50 700 500],'Color', 'w');
    
    figure;
    plot(f_v/1e6,10*log10(psd_v/1e-3),'-r','Linewidth',2); hold on;
    psd_ideal_v = abs(fft(x_ideal_v,NFFT)/n_samples).^2;
    plot(f_v/1e6,10*log10(psd_ideal_v/1e-3),'--b','Linewidth',1.5)
    title('PSD','Interpreter','latex','FontSize', fz);
    legend({'PSD','Ideal'}, 'Interpreter','latex','FontSize', fz-2);
    xlabel('Frequency [MHz]', 'Interpreter','latex','FontSize', fz);
    ylabel('Amplitude [dBm]', 'Interpreter','latex','FontSize', fz);
    grid on; ylim([Ps_dbm-50,Ps_dbm+5])
    set(gcf, 'Position', [50 50 700 500],'Color', 'w');
    
end