clear; close all; clc;

% General
        
config_s.fs_ch = 400e6;
config_s.fs_dsp = 100e6;

config_s.n_fires = 1e3;

config_s.en_noise = 1;
config_s.en_plots = 0;

% TX
config_s.chirp_bw = 100e6;
config_s.chirp_T = 12e-6;
config_s.chirp_P = 100;

config_s.f0 = 24e9;
config_s.pw_tx_dbm = 13;

% Target and channel

config_s.range = 300;
config_s.speed = 70;

config_s.range_max = 300;
config_s.speed_max = 100;

% RX

config_s.snr_db = 10;
config_s.fft_zp = 8;

config_s.n_thr = 20;

%%

figure;
idx_leg = 1;
legends_c = {};
    
for snr_db = 7:2:15
    
    config_s.snr_db = snr_db;

    odata = fmcw_radar_simulator(config_s);
    
    % Theoretical ROC
    [theo_pd_v, theo_pfa_v] = myrocsnr(snr_db);
    p_th = semilogx(theo_pfa_v, theo_pd_v, '--', 'LineWidth', 1.2);
    legends_c{idx_leg} = sprintf("Theo ROC. SNR = %.1f [dB]",snr_db);
    idx_leg = idx_leg + 1;
    hold on;
    
    % Simulated ROC
    p_sim = semilogx(odata.pfa_est_v, odata.pd_est_v, 'o');
    p_sim.Color = p_th.Color;
    p_sim.MarkerFaceColor = p_th.Color;
    p_sim.MarkerEdgeColor = 'k';
    legends_c{idx_leg} = sprintf("Sim ROC. SNR = %.1f [dB]",snr_db);
    idx_leg = idx_leg + 1;
    
end

fz = 16;

% Title, legends, etc
title('FMCW ROC Curves', 'Interpreter','latex','FontSize', fz);

legend(legends_c,'Location','southeast','Interpreter','latex','FontSize', fz-2);
xlabel('Probability of False Alarm', 'Interpreter','latex','FontSize', fz);
ylabel('Probability of Detection', 'Interpreter','latex','FontSize', fz);
grid on
ylim([0,1.1])
xlim([1e-5,1])

set(gcf, 'Position', [50 50 700 700],'Color', 'w');