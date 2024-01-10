clear
close all 
clc

% -- General --
config_s.en_plots = 0;

% -- Tx --
config_s.tx_s.BR = 32e9;                    % Baud rate
config_s.tx_s.M = 16;                       % Cantidad de niveles de la modulacion
config_s.tx_s.NOS = 1;                      % Tasa de sobremuestreo
config_s.tx_s.Lsymbs = 1e4;                 % Cantidad de simbolos
config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
config_s.tx_s.pulse_shaping_type = 0;       % 0: RRC, 1: RC
%---ch---
config_s.ch_awgn.EbNo_db =10; 
config_s.ch_awgn.ISE = 0; %1 activada, 0 desacticada
config_s.ch_awgn.firorder = 17;
config_s.ch_awgn.fc = 20e9;
%ch Configracion de portadora 
config_s.ch_awgn.delta_freq = 50e6; % Offset del LO
config_s.ch_awgn.LW = 00e3; % Ancho de linea [Hz]
config_s.ch_awgn.theta0 = 15/180*pi; % Fase inicial
config_s.ch_awgn.frequency_fluctuations_amp = 0e6;
config_s.ch_awgn.frequency_fluctuations_freq = 1e3;
config_s.ch_awgn.phase_tone_freq = 250e6;
%--PLL--
config_s.rx_s.Kp =0.05;
config_s.rx_s.Ki = config_s.rx_s.Kp/1000;


frec =logspace(log10(1e6), log10(config_s.tx_s.BR/2), 1024);
L=length(frec);
pll_out=zeros(L,1);
parfor idx = 1:length(frec)
    cfg_s = config_s;
    cfg_s.ch_awgn.delta_freq =frec(idx);
    out_c = m_simulatortp6(cfg_s);
    pll_out(idx)= 10*log10(var(out_c.nco_output)/var(out_c.tita_in));
end
position_and_size = [100 100 800 600];
n_freq_pos = 2^16;  
n_freq_frac = 2^16;
n_v = [0, logspace(-log10(n_freq_frac),log10(n_freq_pos),n_freq_pos+n_freq_frac-1)];
f_v = n_v * config_s.tx_s.BR/2/n_freq_pos;
wd_v = n_v * pi/n_freq_pos;

G_th_v = (config_s.rx_s.Kp+config_s.rx_s.Ki.*1./(1-exp(-1j*wd_v))) * 1./(1-exp(-1j.*wd_v));
H_th_v = G_th_v ./ (1+G_th_v.*exp(-1j*wd_v));
H_th_db_v = 20*log10(abs(H_th_v));

figure
semilogx(frec, pll_out,'-r', 'LineWidth', 2);
hold all 
semilogx(f_v, H_th_db_v, '--b', 'LineWidth', 2);
grid on 
xlabel('Frecuencia [Hz]');
ylabel('Amplitude [dB]');
xlim([2e6 4e8]);
title('Jitter Transfer Plot');
xlabel('Time [s]')
ylabel('Amplitudes [dB]')
legend('Simulado', 'Teórico'); % Agrega la leyenda
set(gcf, 'Position',position_and_size,'Color', 'w');



 

