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

%% General configuration

% -- General --
config_s.en_plots = 0;                      % Plots -> ON: 1; OFF: 0

% -- Tx --
config_s.tx_s.BR = 32e9;                    % Baud rate
%config_s.tx_s.M = 16;                      % Cantidad de niveles de la modulacion
config_s.tx_s.NOS = 2;                      % Tasa de sobremuestreo
config_s.tx_s.Lsymbs = 1e5;                 % Cantidad de simbolos
config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
config_s.tx_s.pulse_shaping_type = 0;       % 0: RRC, 1: RC

%% Sweep

M_v = [4 16 64]; %QPSK, QAM16, QAM64
n_M = length(M_v);

out_c = cell(n_M, 1); % Celda para almacenar variables de diferente tipo 

%% Instantiation

for idx_M = 1:n_M   %Barrido para tres indices de modulación
    
    fprintf('- Running %d/%d ...\n', idx_M,n_M)

    cfg_s = config_s;   % Almaceno las config
    cfg_s.tx_s.M = M_v(idx_M); % Agrega el indice M
    
    out_c{idx_M} = m_simulator(cfg_s); % Se envia al main las config

end

file = [out_dir, 'o_data.mat']; % Salve data para posterior proceso
save(file, 'out_c','M_v','config_s');


