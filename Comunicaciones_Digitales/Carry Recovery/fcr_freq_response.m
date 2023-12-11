% Carrier recovery frequency response

clear; close all; clc;

% Baud Rate
BR = 240e9; 

% PLL Kp
kp_v = 5*[1e-4, 2.5e-4, 5e-4, 7.5e-4, 1e-3];

% PLL Ki factor (Factor: Kp/Ki)
ki_factor = 1e4;

% PLL Latency
lat = 128;

% Other variables
fs = BR;
n_freq_pos = 2^16;
n_freq_frac = 2^16;
n_v = [0, logspace(-log10(n_freq_frac),log10(n_freq_pos),n_freq_pos+n_freq_frac-1)];
f_v = n_v * BR/2/n_freq_pos;
wd_v = n_v * pi/n_freq_pos;

% Plot marks
bw_area = [80e6, 100e6];
peak_mark_db = 0.5;

% Plots
fz = 15;

color_c = { [0 0.45 0.74]
            [0.85 0.33 0.1]
            [0.93 0.69 0.13]
            [0.49 0.18 0.56]
            [0.47 0.67 0.19]
            [0.30 0.75 0.93]
            [0.64 0.08 0.18]};

figure

for idx = 1:length(kp_v)

    kp = kp_v(idx);
    ki = kp/ki_factor;

    G_th_v = (kp+ki.*1./(1-exp(-1j*wd_v))) * 1./(1-exp(-1j.*wd_v));
    H_th_v = G_th_v ./ (1+G_th_v.*exp(-1j*lat*wd_v));
    H_th_db_v = 20*log10(abs(H_th_v));

    p = semilogx(f_v, H_th_db_v, '-','Linewidth',1.5);
    p.Color = color_c{idx};
    hold on
    
    leg{idx} = sprintf('Kp = %.2e',kp);

end

plot(xlim, [-3 -3],'-.m','Linewidth',1);
leg{idx+1} = sprintf('-3 dB');
plot(xlim, [peak_mark_db peak_mark_db],'-.c','Linewidth',1);
leg{idx+2} = sprintf('%.2f dB',peak_mark_db);
y_lim = ylim;
a = area(bw_area,[y_lim(1) y_lim(1)]);
a.LineStyle = 'none';
a.FaceAlpha = 0.4;

grid on; ylim([-10, 4]); xlim([1e6, 1e9])
legend(leg,'Location','eo','Interpreter','latex','FontSize', fz-2);
xlabel('Freq [Hz]','Interpreter','latex','FontSize', fz);
ylabel('Amplitude [dB]','Interpreter','latex','FontSize', fz);
tit = sprintf('FCR Jitter Transfer Function. \n BR = %.2f GBd. Ki = Kp/%d. Lat = %d',...
                                                                    BR/1e9,ki_factor, lat);
title(tit,'Interpreter','latex','FontSize', fz);
set(gcf, 'Position', [50 50 700 500],'Color', 'w');