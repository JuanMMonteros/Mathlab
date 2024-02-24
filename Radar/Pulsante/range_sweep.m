clear; close all; clc;

config_s.tau = 5e-9;        % Ancho del pulso [s]
config_s.Po = 7.5e3;        % Potencia instantanea [W]
config_s.max_range = 3e3;   % Rango maximo [m]
config_s.No = 1*(.6e-9)^2;  % PSD del ruido one-side [W/Hz]
config_s.Niters = 1e3;      % Numero de iteraciones (experimentos)

figure;
idx_leg = 1;
legends_c = {};
    
for range = [1800, 2000, 2200, 2500]
    
    config_s.range = range;

    odata = pulsed_radar_simulator(config_s);
    
    % Theoretical ROC
    snr_db = 10*log10(odata.snr_teo);
    [theo_pd_v, theo_pfa_v] = myrocsnr(snr_db);
    p_th = semilogx(theo_pfa_v, theo_pd_v, '--', 'LineWidth', 1.2);
    legends_c{idx_leg} = sprintf("Theo ROC. SNR = %.1f [dB]",snr_db);
    idx_leg = idx_leg + 1;
    hold on;
    
    % Simulated ROC
    p_sim = semilogx(odata.pfa_vector, odata.pd_vector, 'o');
    p_sim.Color = p_th.Color;
    p_sim.MarkerFaceColor = p_th.Color;
    p_sim.MarkerEdgeColor = 'k';
    legends_c{idx_leg} = sprintf("Sim ROC. SNR = %.1f [dB]",snr_db);
    idx_leg = idx_leg + 1;
    
end

fz = 16;

% Title, legends, etc
title('ROC Curves', 'Interpreter','latex','FontSize', fz);

legend(legends_c,'Location','southeast','Interpreter','latex','FontSize', fz-2);
xlabel('Probability of False Alarm', 'Interpreter','latex','FontSize', fz);
ylabel('Probability of Detection', 'Interpreter','latex','FontSize', fz);
grid on
ylim([0,1.1])
xlim([1e-5,1])

set(gcf, 'Position', [50 50 700 700],'Color', 'w');