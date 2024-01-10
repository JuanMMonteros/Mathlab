%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                Test Unitario para el transmisor                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
close all; 
clear
%%
tx_config.M = 16; % Cantidad de niveles de la modulacion QSPK, QAM16, QAM64
tx_config.NOS = 2; % Tasa de sobremuestreo
tx_config.Lsymbs = 100e3; % Cantidad de simbolos
tx_config.rolloff = 0.5; % Rolloff del filtro conformador
tx_config.pulse_shaping_ntaps = 201;    % Num de taps
tx_config.pulse_shaping_type = 0; %Tipo de filtro 0: RRC, 1: RC

odata = TransmisorMQAM(tx_config);  % Se pasa las configuraciones al Tx

tx_symbs = odata.ak_v;  % Obtenemos los simbolos  
tx_data = odata.signal_v; % Datos del transmisor 

%-- Plots --%
scatterplot(tx_symbs)   % Diagrama de constelación
eyediagram(tx_data(1:10e3), tx_config.NOS ) % Diagrama de ojo
pwelch(tx_data) % PSD del Tx
grid on