clear; close all;

%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'out/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

file = [out, 'o_data.mat'];
load(file);

n_fc = length(fc_v);
n_EbNo = length(EbNo_v);

%% Plots

fz = 15;

% BER vs M

ber_theo_v = zeros(n_EbNo,n_fc);
ber_sim_v = zeros(n_EbNo,n_fc);

for idx_fc = 1:n_fc
    
    Fc =  fc_v(idx_fc);
    for idx_EbNo= 1:n_EbNo
    ber_theo_v(idx_EbNo,idx_fc) = out_c{idx_EbNo,idx_fc}.ber_theo;
   ber_sim_v(idx_EbNo,idx_fc)= out_c{idx_EbNo,idx_fc}.ber_sim;
    end
end
%%
figure;
semilogy(EbNo_v, ber_sim_v(:,1), 'r-', 'DisplayName', ['bersim, FC(GHz) = ', num2str(fc_v(1)/1e9)], 'LineWidth', 2);
hold on;
semilogy(EbNo_v, ber_sim_v(:,2), '--r', 'DisplayName', ['bersim, FC(GHz) = ', num2str(fc_v(2)/1e9)], 'LineWidth', 2);
semilogy(EbNo_v, ber_sim_v(:,3), 'b-', 'DisplayName', ['bersim, FC(GHz) = ', num2str(fc_v(3)/1e9)], 'LineWidth', 2);
semilogy(EbNo_v, ber_sim_v(:,4), '--b', 'DisplayName', ['bersim, FC(GHz) = ', num2str(fc_v(4)/1e9)], 'LineWidth', 2);


% for idx_EbNo = 1:n_EbNo
%     plot(M_v, ber_sim_v(idx_EbNo,:), '--o', 'DisplayName', ['bersim, EbNo = ', num2str(EbNo_v(n_EbNo))], 'LineWidth', 1);
% end

tit = sprintf('BER vs EbNo EJ2. BR=%.2fGBd', config_s.tx_s.BR/1e9);
title(tit, 'Interpreter','latex','FontSize', fz);
xlabel('EbNo', 'Interpreter','latex','FontSize', fz);
ylabel('BER', 'Interpreter','latex','FontSize', fz);
ylim([1e-4,0.1]);
legend({},'Location','no','Interpreter','latex','FontSize', fz-2);
set(gcf, 'Position', [50 50 500 500],'Color', 'w');
saveas(gcf,[out_dir,sprintf('figure.png')]);
grid on;
hold on;






% figure;
% p = plot(M_v, ber_theo_v(1), '-o', 'Linewidth', 1);
% hold on;
% grid on;
% plot(EbNo_v, ber_sim_v(1), '--k', 'Linewidth', 1);
% tit = sprintf('BER vs M. BR=%.2fGBd', config_s.tx_s.BR/1e9);
% title(tit, 'Interpreter','latex','FontSize', fz);
% xlabel('M', 'Interpreter','latex','FontSize', fz);
% ylabel('BER', 'Interpreter','latex','FontSize', fz);
% %legend({},'Location','no','Interpreter','latex','FontSize', fz-2);
% set(gcf, 'Position', [50 50 500 500],'Color', 'w');
% %saveas(gcf,[out_dir,sprintf('figure.png')]);