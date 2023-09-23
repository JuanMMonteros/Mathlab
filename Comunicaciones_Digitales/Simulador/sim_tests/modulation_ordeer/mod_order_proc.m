clear; close all;

%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'out/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

file = [out_dir, 'o_data.mat'];
load(file);

n_M = length(M_v);
n_EbNo = length(EbNo_v);

%% Plots

fz = 15;

% BER vs M

ber_theo_v = zeros(n_EbNo,n_M);
ber_sim_v = zeros(n_EbNo,n_M);

for idx_M = 1:n_M
    
    M =  M_v(idx_M);
    for idx_EbNo= 1:n_EbNo
    ber_theo_v(idx_EbNo,idx_M) = out_c{idx_EbNo,idx_M}.ber_theo;
   ber_sim_v(idx_EbNo,idx_M)= out_c{idx_EbNo,idx_M}.ber_sim;
    end
end
%%
figure;
hold on;
for idx_M = 1:n_M
    plot(EbNo_v, ber_theo_v(:,idx_M), 'o-', 'DisplayName', ['bertheo, M = ', num2str(M_v(idx_M))], 'LineWidth', 2);
    plot(EbNo_v, ber_sim_v(:,idx_M), '--o', 'DisplayName', ['bertheo, M = ', num2str(M_v(idx_M))], 'LineWidth', 2);
end
% for idx_EbNo = 1:n_EbNo
%     plot(M_v, ber_sim_v(idx_EbNo,:), '--o', 'DisplayName', ['bersim, EbNo = ', num2str(EbNo_v(n_EbNo))], 'LineWidth', 1);
% end


tit = sprintf('BER vs EbNo. BR=%.2fGBd', config_s.tx_s.BR/1e9);
title(tit, 'Interpreter','latex','FontSize', fz);
xlabel('EbNo', 'Interpreter','latex','FontSize', fz);
ylabel('BER', 'Interpreter','latex','FontSize', fz);
legend({},'Location','no','Interpreter','latex','FontSize', fz-2);
set(gcf, 'Position', [50 50 500 500],'Color', 'w');
saveas(gcf,[out_dir,sprintf('figure.png')]);
grid on;
hold off;






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