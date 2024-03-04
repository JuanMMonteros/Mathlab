clear; close all; clc;

config_s.tau = 5e-9;        % Ancho del pulso [s]
config_s.Po = 7.5e3;        % Potencia instantanea [W]
config_s.max_range = 3e3;   % Rango maximo [m]
config_s.No = 1*(.6e-9)^2;  % PSD del ruido one-side [W/Hz]
config_s.Niters = 100;      % Numero de iteraciones (experimentos)
config_s.NOS = 16;

snr_teo_db = [];
snr_est_db = [];

range_v = [1800, 2000, 2200, 2500];

for range = range_v
    
    config_s.range = range;

    odata = pulsed_radar_simulator(config_s);
    
    % Theoretical ROC
    snr_teo_db(end+1) = 10*log10(odata.snr_teo);
    snr_est_db(end+1) = 10*log10(odata.snr_est);
    
end

figure;
idx_leg = 1;
legends_c = {};

fz = 16;

plot(range_v, snr_teo_db, '--', 'LineWidth', 1.2);
legends_c{idx_leg} = sprintf("Teo");
idx_leg = idx_leg + 1;
hold on;
plot(range_v, snr_est_db, '-xr', 'LineWidth', 1.2);
legends_c{idx_leg} = sprintf("Sim");

% Title, legends, etc
title('SNR vs Range', 'Interpreter','latex','FontSize', fz);

legend(legends_c,'Location','southeast','Interpreter','latex','FontSize', fz-2);
xlabel('Rango [m]', 'Interpreter','latex','FontSize', fz);
ylabel('SNR [dB]', 'Interpreter','latex','FontSize', fz);
grid on

set(gcf, 'Position', [50 50 700 700],'Color', 'w');