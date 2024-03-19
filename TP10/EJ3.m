clear; close all; clc;

config_s.fs_ch =800e6;
config_s.fs_dsp=200e6;

config_s.n_fires=100;

config_s.en_noise=1;
config_s.en_plots=0;

% TX
config_s.chirp_bw=100e6;
config_s.chirp_T=12e-6;
config_s.chirp_P=100;

config_s.f0=23e9;
config_s.pw_tx_dbm=13;

% Target and channel

config_s.range=200;
config_s.speed=70;

config_s.range_max=300;
config_s.speed_max=100;

% RX

config_s.snr_db=10;
config_s.fft_zp=16;

config_s.n_thr=20;


snr=[10,12,14,16,18,20];
n_snr=length(snr);
bw=[50e6,100e6,150e6];
n_bw=length(bw);
T=[6e-6,12e-6,18e-6];
n_T=length(T);
sim_prec_range= zeros(n_snr, n_bw, n_T);
sim_prec_speed= zeros(n_snr, n_bw, n_T);
theo_prec_range=zeros(n_snr, n_bw, n_T);
theo_prec_speed= zeros(n_snr, n_bw, n_T);

for idx_T = 1:n_T 
    for idx_bw = 1:n_bw
        for idx_snr=1:n_snr
            cfg_s=config_s;
            cfg_s.chirp_T=T(idx_T);
            cfg_s.chirp_bw=bw(idx_bw);
            cfg_s.snr_db=snr(idx_snr);
            odata=fmcw_radar_simulator(cfg_s);
            sim_prec_range(idx_snr,idx_bw,idx_T)=odata.range_sim_prec;
            sim_prec_speed(idx_snr,idx_bw,idx_T)=odata.speed_sim_prec;
            theo_prec_range(idx_snr,idx_bw,idx_T)=odata.range_theo_prec;
            theo_prec_speed(idx_snr,idx_bw,idx_T)=odata.speed_theo_prec;
        end
    end
end

%%
%%Plot
for idx_T = 1:n_T
    % Gráfico de precisión del rango vs SNR con escala logarítmica en el eje y
    figure;
    for idx_bw = 1:n_bw
        % Precisión simulada
        semilogy(snr, sim_prec_range(:, idx_bw, idx_T), 'o-', 'LineWidth', 2, 'DisplayName', sprintf('Sim - BW: %d MHz, T: %d µs', bw(idx_bw)/1e6, T(idx_T)*1e6));
        hold on;
        % Precisión teórica
        semilogy(snr, theo_prec_range(:, idx_bw, idx_T), '--', 'LineWidth', 2, 'DisplayName', sprintf('Theo - BW: %d MHz, T: %d µs', bw(idx_bw)/1e6, T(idx_T)*1e6));
    end
    xlabel('SNR (dB)');
    ylabel('Precisión de Rango (m)');
    title(sprintf('Precisión de Rango vs SNR para T = %d µs', T(idx_T)*1e6));
    legend('Location', 'best');
    grid on;

    % Gráfico de precisión de velocidad vs SNR con escala logarítmica en el eje y
    figure;
    for idx_bw = 1:n_bw
        % Precisión simulada
       semilogy(snr, sim_prec_speed(:, idx_bw, idx_T), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Sim - BW: %d MHz, T: %d µs', bw(idx_bw)/1e6, T(idx_T)*1e6));
       hold on;
        % Precisión teórica
        semilogy(snr, theo_prec_speed(:, idx_bw, idx_T), 'o--', 'LineWidth', 2, 'DisplayName', sprintf('Theo - BW: %d MHz, T: %d µs', bw(idx_bw)/1e6, T(idx_T)*1e6));
    end
    xlabel('SNR (dB)');
    ylabel('Precisión de Velocidad (m/s)');
    title(sprintf('Precisión de Velocidad vs SNR para T = %d µs', T(idx_T)*1e6));
    legend('Location', 'best');
    grid on;
end
%%
 file_name = 'out_data.mat';
 save(file_name, 'sim_prec_range','theo_prec_range','sim_prec_speed','theo_prec_speed','config_s','snr','bw','T');