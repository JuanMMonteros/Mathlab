%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% This script shows how to measure the precision of a tone 

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

fft_zp = 32;        % FFT zero pading factor

n_exp = 500;        % Number of experiments

fz = 15;            % Plots Font Size

% SNR vector
snr_db_v = 8:2:24;

%% Processing

n_samples = round(t_meas*fs);

% Time vector
t_v = (0:n_samples-1)*ts;

% dBm to W
Ps = 1e-3*10^(Ps_dbm/10);  % Signal power [W]

% FFT
NFFT = n_samples * fft_zp;
f_v = (0:NFFT-1)*fs/NFFT;

max_fft_v = zeros(1,n_exp);
std_v = zeros(1,length(snr_db_v));
std_theo_v = zeros(1,length(snr_db_v));

for idx_snr = 1:length(snr_db_v)
    
    snr_db = snr_db_v(idx_snr);        % Signal SNR [dB]
    
    snr = 10^(snr_db/10);

    % Noise power: SNR = Es/En = (Ps*t_meas)/(pn/BW);
    % Es: Signal energy = Ps * t_meas
    % En: Noise energy = Pn/BW
    % BW = signal bandwidth. For complex signals BW=fs
    % We can compute the signal power as: Pn = SNR/(Ps*t_meas)*BW

    Pn = (Ps*t_meas)/snr*fs;

    for idx = 1:n_exp

        % Print every 50 experiments 
        if mod(idx,100) == 0
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
        
        [~, idx_max] = max(abs(X_v));
        
        max_fft_v(idx) = (idx_max-1)*fs/NFFT;

    end
    
    std_v(idx_snr) = std(max_fft_v);
    std_theo_v(idx_snr) = 1/(2*pi*t_meas)*sqrt(6/snr);
    
end
   
%% Plots


if en_plots

    figure;
    semilogy(snr_db_v, std_theo_v,'--k','Linewidth',2); hold on;
    semilogy(snr_db_v, std_v,'-xr','Linewidth',2); 
    title('Frequency precision vs SNR','Interpreter','latex','FontSize', fz);
    legend({'CR bound','Simulation'}, 'Interpreter','latex','FontSize', fz-2);
    xlabel('SNR [dB]', 'Interpreter','latex','FontSize', fz);
    ylabel('Frequency precision [Hz]', 'Interpreter','latex','FontSize', fz);
    grid on; 
    set(gcf, 'Position', [50 50 700 500],'Color', 'w');
    
end