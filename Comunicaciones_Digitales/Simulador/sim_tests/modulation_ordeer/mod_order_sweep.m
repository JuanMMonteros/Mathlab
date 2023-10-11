clear

%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'out/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

%% General configuration

 % -- General --
    config_s.en_plots = 0;

    % -- Tx --
    config_s.tx_s.BR = 32e9;                    % Baud rate
    config_s.tx_s.M = 4;                       % Cantidad de niveles de la modulacion
    config_s.tx_s.NOS = 4;                      % Tasa de sobremuestreo
    config_s.tx_s.Lsymbs = 1e6;                 % Cantidad de simbolos
    config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
    config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
    config_s.tx_s.pulse_shaping_type = 0;       % 0: RRC, 1: RC
    %ch
    config_s.ch_awgn.EbNo_db = 3; 
    config_s.ch_awgn.ISE = 1; %1 activada, 0 desacticada
    config_s.ch_awgn.firorder = 17;
    config_s.ch_awgn.fc = 16e9;
    %AGC
    config_s.agc.target = 0.3;
    config_s.agc.adc_phase = 1;
    %Ecualizador Adaptivo
    config_s.ec_s.ntap = 63; 
    config_s.ec_s.N =2;        
    config_s.ec_s.step_cma=2^-11;
    config_s.ec_s.step_dd=2^-11;
    config_s.ec_s.tap_leak_gain=1e-4;
    config_s.ec_s.force_cma_enable=0;
%% Sweep

fc_v = [10 15 20 24];
n_fc = length(fc_v);
EbNo_v=[5 8 10 12];
n_EbNo = length(EbNo_v);

out_c = cell(n_EbNo, n_fc); 

%% Instantiation

for idx_fc = 1:n_fc
    
    fprintf('- Running %d/%d ...\n', idx_fc,n_fc)
    
    for idx_EbNo = 1:n_EbNo
    cfg_s = config_s;
    cfg_s.ch_awgn.EbNo_db =EbNo_v(idx_EbNo);
    cfg_s.ch_awgn.fc =fc_v(idx_fc);
    out_c{idx_EbNo,idx_fc} = m_simulatortp5(cfg_s);
    end

end
%%

file = [out_dir, 'o_data.mat'];
save(file, 'out_c','fc_v','EbNo_v','config_s');
