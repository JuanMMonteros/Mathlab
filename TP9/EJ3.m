clear;
tau_v = [2.5e-9,5e-9,7.5e-9];
range_m(1,:)=1300:100:2200;
range_m(2,:) = range_m(1,:) + 2^(1/4);
range_m(3,:) = range_m(2,:) + 2^(1/4);

for idx_tau = 1:length(tau_v)
    config_s.tau = tau_v(idx_tau);        % Ancho del pulso [s]
    config_s.Po = 5e3;          % Potencia instantanea [W]
    config_s.max_range = 2500;  % Rango maximo [m]
    config_s.No = 1*(.4e-9)^2;  % PSD del ruido one-side [W/Hz]
    config_s.Niters = 1e3;      % Numero de iteraciones (experimentos)
    config_s.NOS=16;
    
    sim_prec = [];
    teo_prec = [];
    snr_db = [];
    
    for range = range_m(idx_tau,:)
        config_s.range=range;
        odata= pulsed_radar_simulator(config_s);
        sim_prec(end+1) = odata.range_sim_prec;
        teo_prec(end+1) = odata.range_theo_prec;
        snr_db(end+1) = 10*log10(odata.snr_teo);
    end
    if idx_tau==1
        figure;
        semilogy(range_m(idx_tau,:),teo_prec,'b--','LineWidth',2);
        hold on;
        grid on;
    end 
 
    semilogy(range_m(idx_tau,:),sim_prec,'-o','LineWidth',2);
end
title(sprintf('Rango VS Precicion con '));
xlabel('Distancia al objetivo (m)');
ylabel('Precicion del Radar (m)');
legend('CR Bound',sprintf("sim tau=%.2e", tau_v(1)),sprintf("sim tau=%.2e", tau_v(2)),sprintf("sim tau=%.2e", tau_v(3)));