clear; close all; clc;

config_s.fs_ch =400e6;
config_s.fs_dsp=100e6;

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
config_s.fft_zp=1;

config_s.n_thr=20;


snr=[8,10,12,14,16,];
n_snr=length(snr);
pfa_th_v= zeros(101,n_snr);
pd_th_v= zeros(101,n_snr);
pd_est_v= zeros(20,n_snr);
pfa_est_v= zeros(20,n_snr);

for idx=1:n_snr
    cfg_s=config_s;
    cfg_s.snr_db = snr(idx);
    odata = fmcw_radar_simulator(cfg_s);
    pfa_th_v(:,idx)= odata.pfa_th_v;
    pd_th_v(:,idx)= odata.pd_th_v;
    pd_est_v(:,idx)= odata.pd_est_v;
    pfa_est_v(:,idx)= odata.pfa_est_v;
end
%%
% -- ROC --
figure;
idx_leg = 1;
legends_c = {};
fz = 15;

for idx = 1:n_snr
    % Theoretical ROC
    semilogx(pfa_th_v(:,idx), pd_th_v(:,idx), '--', 'LineWidth', 1.2);
    legends_c{idx_leg} = sprintf('Theo - SNR %.1f dB', snr(idx));
    idx_leg = idx_leg + 1;
    hold on;

    % Simulated ROC
    semilogx(pfa_est_v(:,idx), pd_est_v(:,idx), 'o-', 'LineWidth', 1.2);
    legends_c{idx_leg} = sprintf('Sim - SNR %.1f dB', snr(idx));
    idx_leg = idx_leg + 1;
end

% Title, legends, etc
tit = "ROC";
title(tit, 'Interpreter','latex','FontSize', fz);

legend(legends_c,'Location','southeast','Interpreter','latex','FontSize', fz-2);
xlabel('Probability of False Alarm', 'Interpreter','latex','FontSize', fz);
ylabel('Probability of Detection', 'Interpreter','latex','FontSize', fz);
grid on
ylim([0,1.1])
xlim([1e-5,1])