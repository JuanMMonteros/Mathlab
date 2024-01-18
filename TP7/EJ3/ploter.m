clear; close all; clc;

%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'out/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

filename = strcat('bsp_16.mat');
file = [out_dir, filename];
load(file);
ebno_ber2e2 = 6.7;
lw=500e3;
file_name = strcat('o_data_', num2str(lw), '.mat');
file = [out_dir, file_name];
load(file);
snr_theo=ebno_ber2e2+ 10*log10(log2(16));
%% semilogxs
fz = 15;
%%
figure
semilogx(filter_lengths, snr_loss - snr_theo , '--o', 'DisplayName', 'BSP 16 fases', 'LineWidth', 2);
hold on;
semilogx(filter_lengths, snr_loss_v(:,1) - snr_theo , 'o-', 'DisplayName', ['payload=', num2str(payload(1))], 'LineWidth', 2);
semilogx(filter_lengths, snr_loss_v(:,2) - snr_theo, 'o-', 'DisplayName', ['payload=', num2str(payload(2))], 'LineWidth', 2);
semilogx(filter_lengths, snr_loss_v(:,3) - snr_theo, 'o-', 'DisplayName', ['payload=', num2str(payload(3))], 'LineWidth', 2);
semilogx(filter_lengths, snr_loss_v(:,4) - snr_theo, 'o-', 'DisplayName', ['payload=', num2str(payload(4))], 'LineWidth', 2);

tit = sprintf('snr_loss vs filter_taps,  LW=%.2fKHz', lw/1e3);
title(tit, 'Interpreter','latex','FontSize', fz);
ylabel('SNR LOSS(db) - BER=2e-2', 'Interpreter','latex','FontSize', fz);
xlabel('filter_taps', 'Interpreter','latex','FontSize', fz);
%ylim([1,-1]);
legend({},'Location','no','Interpreter','latex','FontSize', fz-2);
set(gcf, 'Position', [50 50 500 500],'Color', 'w');
grid on;
hold off;