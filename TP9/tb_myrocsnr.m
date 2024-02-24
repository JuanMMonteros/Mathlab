clear; clc;

snr_db = ;
[pd0,pfa0] = rocsnr(snr_db,'SignalType','NonfluctuatingNoncoherent');
%[pd1,pfa1] = myrocsnr(snr);
[pd1,pfa1] = myrocsnr(snr_db,'MinPfa',1e-10,'MaxPfa',1,'NumPoints',100);

figure; idx_leg = 1; fz = 15;

% ROC 0
semilogx(pfa0, pd0, '-k', 'LineWidth', 2)
legends_c{idx_leg} = "rocsnr()";
idx_leg = idx_leg + 1;
hold on

% ROC 1
p = semilogx(pfa1, pd1, '--r', 'LineWidth', 2);
legends_c{idx_leg} = "myrocsnr()";
idx_leg = idx_leg + 1;

% Title, legends, etc
tit = sprintf('ROC Curves. SNR=%.2f [dB]',snr_db);
title(tit, 'Interpreter','latex','FontSize', fz);

if snr_db < 9
    leg_loc = 'northwest';
else
    leg_loc = 'southeast';
end 

legend(legends_c,'Location',leg_loc,'Interpreter','latex','FontSize', fz-2);
xlabel('Probability of False Alarm', 'Interpreter','latex','FontSize', fz);
ylabel('Probability of Detection', 'Interpreter','latex','FontSize', fz);
grid on
ylim([0,1.1])
xlim([1e-10,1])

set(gcf, 'Position', [50 50 500 500],'Color', 'w');
