%-------------------------------------------------------------------------%
%                                 FULGOR
% Programmer(s): Daniel A. Juarez
% Created on   : October 2023
% Description  : Modulation order and SNR process
%-------------------------------------------------------------------------%

clear 
close all
clc
%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'resultados/'];

figur_dir = mfilename('fullpath');
figur_dir = figur_dir(1:end-length(figur_dir));
figur_dir = [figur_dir, out_dir, 'Figure/']; 

if ~exist(out_dir,'dir')
    mkdir(out_dir);         % Create directory
end
if ~exist(figur_dir,'dir')
    mkdir(figur_dir);
end
%% Load archive .mat 

file = [out_dir, 'cfg.mat'];
load(file);

n_bwCh = length(bwCh_v);
n_ber = length(theo_ber_v);

ebno = config_s.ch_s.EbNo_db; 
rolloff = config_s.tx_s.rolloff;
bw = config_s.tx_s.BR;

ber_est1_m = zeros(n_bwCh,n_ber);
ber_est2_m = zeros(n_bwCh,n_ber);

ber_theo_m = zeros(n_bwCh,n_ber);
ebno_db_m = zeros(n_bwCh,n_ber);

fz = 15;                        % fonts size
color_c = { [0 0.45 0.74]       % Color matrix 
            [0.85 0.33 0.1]
            [0.93 0.69 0.13]
            [0.49 0.18 0.56]
            [0.47 0.67 0.19]
            [0.30 0.75 0.93]
            [0.64 0.08 0.18]};


%% -- BER vs EbNo Opcion: filter(signal) + noise -- %
opc = 0;
path = mfilename('fullpath');
path = path(1:end-length(path)); 

for idx_1 = 1:n_bwCh
    
    bwCh = bwCh_v(idx_1);
    archive_mat = sprintf('FC%.2e',bwCh);
    name = sprintf('BW%.2e[Hz]_EbNo%d[dB]_FC%.2e[Hz]_Opc%d', bw, ebno, ...
        bwCh, opc);
    
    path = [out_dir, name, '/'];    % Load out_Mx.mat 
    file = [path, 'out_',archive_mat,'.mat'];
    load(file);

    for idx_2 = 1:n_ber

        ber_est1_m(idx_1,idx_2) = out_c{idx_2}.ber_est;
        ber_theo_m(idx_1,idx_2) = out_c{idx_2}.ber_theo;

    end
end

%% -- Plot --
figure; 
leg = {}; 
idx_leg = 1;
        
for idx_1 = 1:n_bwCh
    
    bwCh = bwCh_v(idx_1);
    ebno_db_v = get_ebno_from_theo_ber(theo_ber_v, config_s.tx_s.M);
    
    ber_theo_v = ber_theo_m(idx_1,:);
    ber_est1_v = ber_est1_m(idx_1,:);
    

    p = semilogy(ebno_db_v, ber_est1_v, '-o', 'Linewidth', 1);
    p.MarkerFaceColor = color_c{idx_1};
    p.MarkerEdgeColor = 'k';
    p.Color = color_c{idx_1}; 
    hold on; grid on;
    leg{idx_leg} = sprintf('BER Est. FC=%.2e[Hz]',bwCh); 
    idx_leg = idx_leg + 1;
        
end

p = semilogy(ebno_db_v, ber_theo_v, '--', 'Linewidth', 1.5);
p.Color = 'k'; 
leg{idx_leg} = sprintf('BER Theo. bw=%.2e[Hz]',bw); 
idx_leg = idx_leg + 1;

xlabel('EbNo [GBd]', 'Interpreter','latex','FontSize', fz);
ylabel('BER', 'Interpreter','latex','FontSize', fz);
legend(leg,'Location','sw','Interpreter','latex','FontSize', fz-2);
grid on; 
ylim([10e-4,10e-1]);

tit = ['BER vs EbNo. ',...
        sprintf('BR=%.0f[GHz]', config_s.tx_s.BR/1e9)];
title(tit, 'Interpreter','latex','FontSize', fz);
set(gcf, 'Position', [1200 500 1100 700],'Color', 'w');
saveas(gcf,[figur_dir,sprintf('BERvsEbNo_EbNo%d_OPC%d.png', ...
            ebno,opc)]);


%% -- BER vs EbNo Opcion: signal + noise -- 
opc = 1;
path = mfilename('fullpath');
path = path(1:end-length(path)); 

for idx_1 = 1:n_bwCh
    
    bwCh = bwCh_v(idx_1);
    archive_mat = sprintf('FC%.2e',bwCh);
    name = sprintf('BW%.2e[Hz]_EbNo%d[dB]_FC%.2e[Hz]_Opc%d', bw, ebno, ...
        bwCh, opc);
    
    path = [out_dir, name, '/'];    % Load out_Mx.mat 
    file = [path, 'out_',archive_mat,'.mat'];
    load(file);

    for idx_2 = 1:n_ber

        ber_est2_m(idx_1,idx_2) = out_c{idx_2}.ber_est;
        ber_theo_m(idx_1,idx_2) = out_c{idx_2}.ber_theo;

    end
end

%% -- Plot --
figure; 
leg = {}; 
idx_leg = 1;

for idx_1 = 1:n_bwCh

    bwCh = bwCh_v(idx_1);
    ebno_db_v = get_ebno_from_theo_ber(theo_ber_v, config_s.tx_s.M);

    ber_theo_v = ber_theo_m(idx_1,:);
    ber_est2_v = ber_est2_m(idx_1,:);


    p = semilogy(ebno_db_v, ber_est2_v, '-o', 'Linewidth', 1);
    p.MarkerFaceColor = color_c{idx_1};
    p.MarkerEdgeColor = 'k';
    p.Color = color_c{idx_1}; 
    hold on; grid on;
    leg{idx_leg} = sprintf('BER Est. FC=%.2e[Hz]',bwCh); 
    idx_leg = idx_leg + 1;

end

p = semilogy(ebno_db_v, ber_theo_v, '--', 'Linewidth', 1.5);
p.Color = 'k'; 
leg{idx_leg} = sprintf('BER Theo. bw=%.2e[Hz]',bw); 
idx_leg = idx_leg + 1;

xlabel('EbNo [GBd]', 'Interpreter','latex','FontSize', fz);
ylabel('BER', 'Interpreter','latex','FontSize', fz);
legend(leg,'Location','sw','Interpreter','latex','FontSize', fz-2);
grid on; 
ylim([10e-4,10e-1]);

tit = ['BER vs EbNo. ',...
        sprintf('BR=%.0f[GHz]', config_s.tx_s.BR/1e9)];
title(tit, 'Interpreter','latex','FontSize', fz);
set(gcf, 'Position', [1200 500 1100 700],'Color', 'w');
saveas(gcf,[figur_dir,sprintf('BERvsEbNo_EbNo%d_OPC%d.png', ...
            ebno,opc)]);


%% -- Plot ambos casos --
figure; 
leg = {}; 
idx_leg = 1;

for idx_opc = 1:2

    % ber_theo_v = ber_theo_m(idx_1,:);

    if idx_opc == 1

        for idx_1 = 1:n_bwCh
            bwCh = bwCh_v(idx_1);
            ebno_db_v = get_ebno_from_theo_ber(theo_ber_v, config_s.tx_s.M);
    
            ber_est1_v = ber_est1_m(idx_1,:);
    
            p = semilogy(ebno_db_v, ber_est1_v, '-o', 'Linewidth', 1);
            p.MarkerFaceColor = color_c{idx_1};
            p.MarkerEdgeColor = 'k';
            p.Color = color_c{idx_1};
            hold on; grid on;
            leg{idx_leg} = sprintf('BER_Est_pre_noise. FC=%.2e[Hz]',bwCh); 
            idx_leg = idx_leg + 1;        
        end
    end
    if idx_opc == 2
        for idx_1 = 1:n_bwCh
            bwCh = bwCh_v(idx_1);
            ebno_db_v = get_ebno_from_theo_ber(theo_ber_v, config_s.tx_s.M);
    
            ber_est2_v = ber_est2_m(idx_1,:);
    
            p = semilogy(ebno_db_v, ber_est2_v, '--', 'Linewidth', 1);
            p.MarkerFaceColor = color_c{idx_1};
            p.MarkerEdgeColor = 'k';
            p.Color = color_c{idx_1};
            hold on; grid on;
            leg{idx_leg} = sprintf('BER_Est_post_noise. FC=%.2e[Hz]',bwCh); 
            idx_leg = idx_leg + 1;        
        end
    end

end

xlabel('EbNo [GBd]', 'Interpreter','latex','FontSize', fz);
ylabel('BER', 'Interpreter','latex','FontSize', fz);
legend(leg,'Location','sw','Interpreter','latex','FontSize', fz-2);
grid on; 
ylim([10e-4,5e-1]);

tit = ['BER vs EbNo. ',...
        sprintf('BR=%.0f[GHz]', config_s.tx_s.BR/1e9)];
title(tit, 'Interpreter','latex','FontSize', fz);
set(gcf, 'Position', [1200 500 1100 700],'Color', 'w');
saveas(gcf,[figur_dir,sprintf('BERvsEbNo_AmbasCurvas.png')]);



