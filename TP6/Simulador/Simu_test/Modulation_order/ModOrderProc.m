%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     SIMULADOR QAM-M COMPLETO                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
close all; 
clear
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

%% Plots

fz = 15;
 
% color_c = { [0 0.45 0.74]
%             [0.85 0.33 0.1]
%             [0.93 0.69 0.13]
%             [0.49 0.18 0.56]
%             [0.47 0.67 0.19]
%             [0.30 0.75 0.93]
%             [0.64 0.08 0.18]};

% BER vs M

ber_v = zeros(n_M,1);

for idx_M = 1:n_M
    
    M =  M_v(idx_M);
    
    ber_v(idx_M) = out_c{idx_M}.ber_theo;
        
end 

figure;
p = plot(M_v, ber_v, '-o', 'Linewidth', 1);
% p.MarkerFaceColor = color_c{1};
% p.MarkerEdgeColor = 'k';
tit = sprintf('BER vs M. BR=%.2fGBd', config_s.tx_s.BR/1e9);
title(tit, 'Interpreter','latex','FontSize', fz);
xlabel('M', 'Interpreter','latex','FontSize', fz);
ylabel('BER', 'Interpreter','latex','FontSize', fz);
%legend({},'Location','no','Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [50 50 500 500],'Color', 'w');
saveas(gcf,[out_dir,sprintf('figure.png')]);


