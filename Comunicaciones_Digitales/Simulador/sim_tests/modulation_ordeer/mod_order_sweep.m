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
    config_s.tx_s.M = 16;                       % Cantidad de niveles de la modulacion
    config_s.tx_s.NOS = 2;                      % Tasa de sobremuestreo
    config_s.tx_s.Lsymbs = 1e6;                 % Cantidad de simbolos
    config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
    config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
    config_s.tx_s.pulse_shaping_type = 0;       % 0: RRC, 1: RC
    config_s.ch_awgn.EbNo_db = 8; 
    config_s.ch_awgn.M =  config_s.tx_s.M;                       % Cantidad de niveles de la modulacion
    config_s.ch_awgn.NOS =  config_s.tx_s.NOS;
    config_s.rx_s.filter_type = 1 ;        % 1: MF, 2: impulso
    config_s.rx_s.NOS =  config_s.tx_s.NOS;
    config_s.rx_s.ntaps = config_s.tx_s.pulse_shaping_ntaps;
    config_s.ber_s.M = config_s.tx_s.M;
    config_s.ber_s.EbNo_db = config_s.ch_awgn.EbNo_db; 

%% Sweep

M_v = [4 16];
n_M = length(M_v);
EbNo_v=[5 8 10 12];
n_EbNo = length(EbNo_v);

out_c = cell(n_EbNo, n_M); 

%% Instantiation

for idx_M = 1:n_M
    
    fprintf('- Running %d/%d ...\n', idx_M,n_M)
    
    for idx_EbNo = 1:n_EbNo
    cfg_s = config_s;
    cfg_s.tx_s.M = M_v(idx_M);
    cfg_s.ch_awgn.M = M_v(idx_M);
    cfg_s.ber_s.M = M_v(idx_M);
    cfg_s.ch_awgn.EbNo_db =EbNo_v(idx_EbNo);
    cfg_s.ber_s.EbNo_db =EbNo_v(idx_EbNo);
    out_c{idx_EbNo,idx_M} = m_simulator(cfg_s);
    end

end
%%

file = [out_dir, 'o_data.mat'];
save(file, 'out_c','M_v','EbNo_v','config_s');
