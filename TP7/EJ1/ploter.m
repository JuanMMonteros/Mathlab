clear; close all; clc;

%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'out/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
lw=[0,10e3,10e4,10e4,50e4];
for idx_lw = 1:length(lw)
file_name = strcat('o_data_', num2str(lw(idx_lw)), '.mat');
file = [out_dir, file_name];
load(file);

%% Plots
fz = 15;
%%
figure
plot(filter_lengths, snr_loss_v(:,1), 'o-', 'DisplayName', ['n\_phases=', num2str(phases(1))], 'LineWidth', 2);
hold on;
plot(filter_lengths, snr_loss_v(:,2), 'o-', 'DisplayName', ['n\_phases=', num2str(phases(2))], 'LineWidth', 2);
plot(filter_lengths, snr_loss_v(:,3), 'o-', 'DisplayName', ['n\_phases=', num2str(phases(3))], 'LineWidth', 2);
plot(filter_lengths, snr_loss_v(:,4), 'o-', 'DisplayName', ['n\_phases=', num2str(phases(4))], 'LineWidth', 2);


% for idx_EbNo = 1:n_EbNo
%     plot(M_v, ber_sim_v(idx_EbNo,:), '--o', 'DisplayName', ['bersim, EbNo = ', num2str(EbNo_v(n_EbNo))], 'LineWidth', 1);
% end

tit = sprintf('snr_loss vs filter_lengths. LW=%.2fKHz', lw(idx_lw)/1e3);
title(tit, 'Interpreter','latex','FontSize', fz);
ylabel('SNR LOSS', 'Interpreter','latex','FontSize', fz);
xlabel('filter_lengths', 'Interpreter','latex','FontSize', fz);
%ylim([1e-4,0.1]);
legend({},'Location','no','Interpreter','latex','FontSize', fz-2);
set(gcf, 'Position', [50 50 500 500],'Color', 'w');
grid on;
hold off;
end