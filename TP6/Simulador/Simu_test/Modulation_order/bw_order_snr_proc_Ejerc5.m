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

file = [out_dir, 'cfg2.mat'];
load(file);

n_taps = length(taps_v);
n_ber = length(theo_ber_v);
n_step_dd = length(step_DD_v);

ebno = config_s.ch_s.EbNo_db; 
rolloff = config_s.tx_s.rolloff;
bw = config_s.tx_s.BR; 

ber_est1_m = zeros(n_taps,n_ber);
ber_est2_m = zeros(n_taps,n_ber);
ber_est3_m = zeros(n_taps,n_ber);
ber_est4_m = zeros(n_taps,n_ber);
ber_theo1_m = zeros(n_taps,n_ber);
ber_theo2_m = zeros(n_taps,n_ber);
ber_theo3_m = zeros(n_taps,n_ber);
ber_theo4_m = zeros(n_taps,n_ber);
ebno_db_m = zeros(n_taps,n_ber);

fz = 15;                        % fonts size
color_c = { [0 0.45 0.74]       % Color matrix 
            [0.85 0.33 0.1]
            [0.93 0.69 0.13]
            [0.49 0.18 0.56]
            [0.47 0.67 0.19]
            [0.30 0.75 0.93]
            [0.64 0.08 0.18]};


%% Read data

path = mfilename('fullpath');
path = path(1:end-length(path)); 

for idx_step = 1:n_step_dd
    step_DD = step_DD_v(idx_step);

    for idx_1 = 1:n_taps
        
        taps = taps_v(idx_1);
        archive_mat = sprintf('Taps%d', taps);
        name = sprintf('BW%.2e[Hz]_Taps%d_Step_DD%.2e',bw,taps,step_DD);
        
        path = [out_dir, name, '/'];    % Load out_Mx.mat 
        file = [path, 'out_',archive_mat,'.mat'];
        load(file);
        
        for idx_2 = 1:n_ber
            if step_DD == 2^-9
                ber_est1_m(idx_1,idx_2) = out_c{idx_2}.ber_est;
                ber_theo1_m(idx_1,idx_2) = out_c{idx_2}.ber_theo;
            end
            if step_DD == 4^-9
                ber_est2_m(idx_1,idx_2) = out_c{idx_2}.ber_est;
                ber_theo2_m(idx_1,idx_2) = out_c{idx_2}.ber_theo;
            end
            if step_DD == 6^-9
                ber_est3_m(idx_1,idx_2) = out_c{idx_2}.ber_est;
                ber_theo3_m(idx_1,idx_2) = out_c{idx_2}.ber_theo;
            end
            if step_DD == 8^-9
                ber_est4_m(idx_1,idx_2) = out_c{idx_2}.ber_est;
                ber_theo4_m(idx_1,idx_2) = out_c{idx_2}.ber_theo;
            end            
        end
    end
end


%%
%  -- BER vs EbNo Opcion: filter(signal) + noise -- %
figure; 
leg = {}; 
idx_leg = 1;
Step_DD = step_DD_v(4);


for idx_1 = 1:n_taps

    ber_theo_v = ber_theo4_m(idx_1,:);
    ber_est_v = ber_est4_m(idx_1,:);

    taps = taps_v(idx_1);
    ebno_db_v = get_ebno_from_theo_ber(theo_ber_v, config_s.tx_s.M);
        
    p = semilogy(ebno_db_v, ber_est_v, '-o', 'Linewidth', 1);
    p.MarkerFaceColor = color_c{idx_1};
    p.MarkerEdgeColor = 'k';
    p.Color = color_c{idx_1};
    hold on; grid on;
    leg{idx_leg} = sprintf('Est. bw=%.2e, Taps%d',bw, taps); 
    idx_leg = idx_leg + 1;

end

p = semilogy(ebno_db_v, ber_theo_v, '--', 'Linewidth', 1.5);
p.Color = 'k'; 
hold on; grid on;
leg{idx_leg} = sprintf('Theo. bw=%.2e',bw); 
idx_leg = idx_leg + 1;

xlabel('Taps', 'Interpreter','latex','FontSize', fz);
ylabel('BER', 'Interpreter','latex','FontSize', fz);
legend(leg,'Location','sw','Interpreter','latex','FontSize', fz-2);
grid on; 
ylim([3e-4,0.2]);

tit = ['BER vs EbNo. ',...
        sprintf('Step DD=%.2e', Step_DD)];
title(tit, 'Interpreter','latex','FontSize', fz);
set(gcf, 'Position', [1200 500 1100 700],'Color', 'w');
saveas(gcf,[figur_dir,sprintf('BERvsEbNo_Step_DD%.2e.png', ...
            Step_DD)]);


%%
% -- SNR Penalty vs Bw Opcion: filter(signal) + noise  -- %
ber_int = 1*10e-3;
snr_loss1_db_v = zeros(n_taps,1);
snr_loss2_db_v = zeros(n_taps,1);
snr_loss3_db_v = zeros(n_taps,1);
snr_loss4_db_v = zeros(n_taps,1);
taps_v = [201 501 801];    % Bw order in [GHz]
n_napts = length(taps_v);
name = {};
idx_name = 1;

for idx_1 = 1:n_taps

    ebno_db_v = get_ebno_from_theo_ber(theo_ber_v,config_s.tx_s.M);
    ber_theo_v = ber_theo1_m(idx_1,:);
    ber_est_v = ber_est1_m(idx_1,:);

    ebno_sim_db = interp1(log10(ber_est_v), ebno_db_v, log10(ber_int));
    ebno_theo_db = interp1(log10(ber_theo_v), ebno_db_v, log10(ber_int));

    snr_loss1_db_v(idx_1) =  ebno_sim_db - ebno_theo_db;
    
end

for idx_1 = 1:n_taps

    ebno_db_v = get_ebno_from_theo_ber(theo_ber_v,config_s.tx_s.M);
    ber_theo_v = ber_theo2_m(idx_1,:);
    ber_est_v = ber_est2_m(idx_1,:);

    ebno_sim_db = interp1(log10(ber_est_v), ebno_db_v, log10(ber_int));
    ebno_theo_db = interp1(log10(ber_theo_v), ebno_db_v, log10(ber_int));

    snr_loss2_db_v(idx_1) =  ebno_sim_db - ebno_theo_db;

end
for idx_1 = 1:n_taps

    ebno_db_v = get_ebno_from_theo_ber(theo_ber_v,config_s.tx_s.M);
    ber_theo_v = ber_theo3_m(idx_1,:);
    ber_est_v = ber_est3_m(idx_1,:);

    ebno_sim_db = interp1(log10(ber_est_v), ebno_db_v, log10(ber_int));
    ebno_theo_db = interp1(log10(ber_theo_v), ebno_db_v, log10(ber_int));

    snr_loss3_db_v(idx_1) =  ebno_sim_db - ebno_theo_db;

end
for idx_1 = 1:n_taps

    ebno_db_v = get_ebno_from_theo_ber(theo_ber_v,config_s.tx_s.M);
    ber_theo_v = ber_theo4_m(idx_1,:);
    ber_est_v = ber_est4_m(idx_1,:);

    ebno_sim_db = interp1(log10(ber_est_v), ebno_db_v, log10(ber_int));
    ebno_theo_db = interp1(log10(ber_theo_v), ebno_db_v, log10(ber_int));

    snr_loss4_db_v(idx_1) =  ebno_sim_db - ebno_theo_db;

end

figure; 

p = plot(taps_v, snr_loss1_db_v, '-o', 'Linewidth', 2);
p.MarkerFaceColor = 'r';
p.MarkerEdgeColor = 'k';
p.Color = 'r'; 
grid on; hold on;

p = plot(taps_v, snr_loss2_db_v, '-o', 'Linewidth', 2);
p.MarkerFaceColor = 'b';
p.MarkerEdgeColor = 'k';
p.Color = 'b'; 

p = plot(taps_v, snr_loss3_db_v, '-o', 'Linewidth', 2);
p.MarkerFaceColor = 'y';
p.MarkerEdgeColor = 'k';
p.Color = 'y'; 

p = plot(taps_v, snr_loss4_db_v, '-o', 'Linewidth', 2);
p.MarkerFaceColor = 'g';
p.MarkerEdgeColor = 'k';
p.Color = 'g'; 

xlabel('Taps LMS', 'Interpreter','latex','FontSize', fz);
ylab = sprintf('SNR loss [dB] @ BER = %.1e',ber_int);
ylabel(ylab, 'Interpreter','latex','FontSize', fz);
 
tit = 'SNR loss vs Taps. ';
title(tit, 'Interpreter','latex','FontSize', fz);

for idx_step = 1:n_step_dd
    name{idx_name} = sprintf('Step DD %.2e',step_DD_v(idx_step));
    idx_name = idx_name+1; 
end    
legend(name,'Location','sw','Interpreter','latex','FontSize',fz-2);

set(gcf, 'Position', [1200 500 1100 700],'Color', 'w');
saveas(gcf,[figur_dir,sprintf('SNRLossVsNtaps.png')]);
