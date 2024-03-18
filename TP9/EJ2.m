clear; close all; clc;

config_s.tau = 5e-9;        % Ancho del pulso [s]
config_s.Po = 5e3;          % Potencia instantanea [W]
config_s.max_range = 2500;  % Rango maximo [m]
config_s.No = 1*(.4e-9)^2;  % PSD del ruido one-side [W/Hz]
config_s.Niters = 1e3;      % Numero de iteraciones (experimentos)
config_s.NOS=16;
range = [1000,1400,1800, 2000, 2200, 2500];
n_range = length(range);
snr_teo = zeros(1, n_range);
snr_sim = zeros(1, n_range);

for idx = 1:n_range
    config_s.range = range(idx);
    odata = pulsed_radar_simulator(config_s);
    snr_teo(idx) = 10*log10(odata.snr_teo);
    snr_sim(idx) = 10*log10(odata.snr_est);
%     snr_teo(idx) = odata.snr_teo;
%     snr_sim(idx) = odata.snr_est;
end


figure;
hold on;
plot(range, snr_sim, 'o-', 'LineWidth', 2);  % snr_sim con una línea rojo continuo
plot(range, snr_teo, 'o--', 'LineWidth', 2); % snr_teo con una línea discontua azul
xlabel('Rango (m)');
ylabel('SNR (dB)');
xlim([1000,2500]);
title('SNR vs Distancia del Target', 'FontSize', 16);
legend('SNR sim', 'SNR teo');
grid on;
hold off;

