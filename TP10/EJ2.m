clear; close all; clc;

config_s.fs_ch =400e6;
config_s.fs_dsp=100e6;

config_s.n_fires=1000;

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
config_s.fft_zp=8;

config_s.n_thr=20;

snr=[12,14,16,18];
n_snr = length(snr);
snr_sim = zeros(1, n_snr);

parfor idx = 1:n_snr
    config_s.snr_db = snr(idx);
    odata = fmcw_radar_simulator(config_s);
    snr_sim(idx) = odata.snr_est_db;
end

figure;
plot(snr, snr_sim, 'o-', 'LineWidth', 2);  
xlabel('SNR Simulada(dB)');
ylabel('SNR Teorica(dB)');
xlim([12,18]);
title('SNR Teorica vs Simulada', 'FontSize', 16);
grid on;

