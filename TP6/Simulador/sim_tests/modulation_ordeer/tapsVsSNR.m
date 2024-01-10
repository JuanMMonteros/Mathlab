clear
close all 
clc
%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'out/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

%% General configuration

    config_s.en_plots = 0;

    % -- Tx --
    config_s.tx_s.BR = 32e9;                    % Baud rate
    config_s.tx_s.M = 4;                       % Cantidad de niveles de la modulacion
    config_s.tx_s.NOS = 2;                      % Tasa de sobremuestreo
    config_s.tx_s.Lsymbs = 1e6;                 % Cantidad de simbolos
    config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
    config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
    config_s.tx_s.pulse_shaping_type = 0;       % 0: RRC, 1: RC
    %ch
    config_s.ch_awgn.EbNo_db =10; 
    config_s.ch_awgn.ISE = 1; %1 activada, 0 desacticada
    config_s.ch_awgn.firorder = 17;
    config_s.ch_awgn.fc = 20e9;
    %AGC
    config_s.agc.target = 0.3;
    config_s.agc.adc_phase = 1;
    %Ecualizador Adaptivo
    config_s.ec_s.ntaps = 63; 
    config_s.ec_s.N =2;        
    config_s.ec_s.step_cma=2^-11;
    config_s.ec_s.step_dd=2^-11;
    config_s.ec_s.tap_leak_gain=1e-7;
    config_s.ec_s.force_cma_enable=0;
    %% Sweep

taps_v = [1 5 31 63];
n_taps = length(taps_v);
EbNo_v=[0 3 4 5 6 8 10 14];
n_EbNo = length(EbNo_v);

out_c = cell(n_EbNo, n_taps); 

%% Instantiation

for idx_taps = 1:n_taps
    
    fprintf('- Running %d/%d ...\n', idx_taps,n_taps)
    
    parfor idx_EbNo = 1:n_EbNo
    cfg_s = config_s;
    cfg_s.ch_awgn.EbNo_db =EbNo_v(idx_EbNo);
    cfg_s.ec_s.ntaps =taps_v(idx_taps);
    out_c{idx_EbNo,idx_taps} = m_simulatortp5(cfg_s);
    fprintf('/////- bertheo %d \n - bersim %d \n', out_c{idx_EbNo,idx_taps}.ber_theo,out_c{idx_EbNo,idx_taps}.ber_sim)
    end

end
%%

file = [out_dir, 'o_data.mat'];
save(file, 'out_c','taps_v','EbNo_v','config_s');
