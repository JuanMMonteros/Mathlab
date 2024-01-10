%-------------------------------------------------------------------------%
%                                 FULGOR
% Programmer(s): Daniel A. Juarez
% Created on   : October 2023
% Description  : Bandwidth and SNR sweep
%-------------------------------------------------------------------------%

clear 
close all
clc
%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'resultados/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

%% General configuration

% -- General --
config_s.en_plots = 0;                      % Plots -> ON: 1; OFF: 0

% -- Tx -- % 
config_s.tx_s.BR = 32e9;                    % Baud rate
config_s.tx_s.M = 4;                        % Niveles de la modulacion
config_s.tx_s.NOS = 2;                      % Tasa de sobremuestreo
config_s.tx_s.Lsymbs = 1e5;                 % Cantidad de simbolos
config_s.tx_s.rolloff = 0.5;                % Rolloff filtro conformador
config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
config_s.tx_s.pulse_shaping_type = 0;       % 0: RRC, 1: RC

% -- Ch --
config_s.ch_s.enable_noise = 1;             % noise -> 0:OFF ; 1:ON
config_s.ch_s.EbNo_db = 3;                  % EbNo [dB]
config_s.ch_s.enable_Bw_limit = 1;          % BW limit -> 0:OFF ; 1:ON
config_s.ch_s.bw_order = 17;                % BW order
config_s.ch_s.bw = 20e9;                    % BW lim [Hz]
config_s.ch_s.opc = 0;                      % 1: signal + noise -> filter
                                            % 0: filter(signal) -> filter + noise                            
% -- Rx FSE -- 
config_s.rx_s.Ntaps_lms = 63;               % Num de taps   
config_s.rx_s.step_cma = 2e-10;             % Tap para el CMA
config_s.rx_s.step_dd = 2^-9;               % Tap para DD   
config_s.rx_s.tap_leak = 10e-3;             % Paso del leakeag
config_s.rx_s.force_cma_enable = 0;         % Force CMA, ON -> 1, OFF -> 0 

%% Sweep

bwCh_v = [16e9 18e9 20e9];    % Bw order in [GHz]
n_bwCh = length(bwCh_v);

ebno = config_s.ch_s.EbNo_db; 
rolloff = config_s.tx_s.rolloff;
bw = config_s.tx_s.BR;

ber_max = 1e-1; % BER limit
ber_min = 1e-4;
n_ber = 10;
linespace_v = linspace(log10(ber_min), log10(ber_max), n_ber);
theo_ber_v = 10.^(linespace_v); % BER theo

% Save configuration in cell
file = [out_dir, 'cfg.mat'];
save(file, 'bwCh_v','theo_ber_v','config_s');

out_c = cell(n_ber, 1); 

%% Barrido para la opción: filter(signal) -> filter + noise 

path = mfilename('fullpath');   
path = out_dir(1:end-length(path));
opcCh = 0;   

for idx_bw = 1:n_bwCh
    
    bwCh = bwCh_v(idx_bw); % Vars 
    ebno_db_v = get_ebno_from_theo_ber(theo_ber_v,config_s.tx_s.M);
    config_s.ch_s.bw = bwCh;  % Se barre el Bw del canal 
    
    for idx_ber = 1:n_ber        
        config_s.ch_s.EbNo_db = ebno_db_v(idx_ber);   
        fprintf('- Running OPC%d %d/%d-%d/%d ...\n',opcCh, idx_bw,n_bwCh,idx_ber,n_ber)
        out_c{idx_ber} = m_simulator(config_s);
    end

    name = sprintf('FC%.2e',bwCh);          % Generate filename 

    % Generate dir name 
    direc = sprintf('BW%.2e[Hz]_EbNo%d[dB]_FC%.2e[Hz]_Opc%d', bw, ebno, bwCh, opcCh); 
    path = [out_dir, direc,'/'];      % Path directory

    if ~exist(path,'dir')
    mkdir(path);                      % Create directory
    end

    file = [path,'out_',name,'.mat']; % save data in direc
    save(file, 'out_c');
    
end

%% Barrido para la opción: signal + noise -> filter

config_s.ch_s.opc = 1;
opcCh  = 1;

path = mfilename('fullpath');   
path = out_dir(1:end-length(path));
   
for idx_bw = 1:n_bwCh
    
    bwCh = bwCh_v(idx_bw); % Vars 
    ebno_db_v = get_ebno_from_theo_ber(theo_ber_v,config_s.tx_s.M);
    config_s.ch_s.bw = bwCh;  % Se barre el Bw del canal 
    
    for idx_ber = 1:n_ber
        config_s.ch_s.EbNo_db = ebno_db_v(idx_ber);   
        fprintf('- Running OPC%d %d/%d-%d/%d ...\n',opcCh, idx_bw,n_bwCh,idx_ber,n_ber)
        out_c{idx_ber} = m_simulator(config_s);
    end

    name = sprintf('FC%.2e',bwCh);    % Generate filename 

    % Generate dir name 
    direc = sprintf('BW%.2e[Hz]_EbNo%d[dB]_FC%.2e[Hz]_Opc%d', bw, ebno, bwCh,opcCh); 
    path = [out_dir, direc,'/'];      % Path directory

    if ~exist(path,'dir')
    mkdir(path);                      % Create directory
    end

    file = [path,'out_',name,'.mat']; % save data in direc
    save(file, 'out_c');
    
end
