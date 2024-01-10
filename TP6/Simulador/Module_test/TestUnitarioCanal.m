%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Test Unitario para el canal                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
close all; 
clear
%% -- SETTINGS --

ch_config.BR = 32e9;            % Niveles de modulación
ch_config.M = 16;               % Orden de modulación
ch_config.NOS = 2;              % Factor de sobremuestreo
ch_config.enable_noise = 0;     % noise -> 0:OFF ; 1:ON
ch_config.EbNo_db = 8;          % EbNo [dB]
ch_config.enable_Bw_limit = 0;  % BW limit -> 0:OFF ; 1:ON
ch_config.bw_order = 8;         % BW order
ch_config.bw = 10e9;            % BW lim [Hz]


%_________Generate signal_____________%

% -- settings TX -- %
M = ch_config.M;                % Order modul
BR = ch_config.BR;              % Baud Rate
N = ch_config.NOS;              % Oversampling rate
EbNo_db = ch_config.EbNo_db;    % EbNo in dB
L = 1e5;                        % Simulation Length
rolloff = 0.5;                  % Pulse shaping rolloff
h_taps = 201;                   % Pulse shaping taps


% -- basic Tx -- %
fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 consecutive samples at Tx output

% generate labels
dec_labels = randi([0 M-1], L, 1); 
% generate simbol
ak = qammod(dec_labels,M); 
% Upsampling to change sampling rate
xup = N * upsample(ak,N);
% Filter to interpolate the signal
h = RootRainCosine(BR/2, fs, rolloff, h_taps, 0);
% output filter tx
yup = filter(h,1,xup);

% -- PROCESS -- %
i_ch_s.signal_v = yup;
odata = CanalDeTransmision(i_ch_s, ch_config); 

signal_v = odata.signal_v; % Signal ouput filter ch
hch_v    = odata.hch_v;   % filter ch

%-- Plots --%
scatterplot(signal_v)
grid on
title('Constelación salida del filtro de canal');
set(gcf, 'Position', [50 50 600 600],'Color', 'w');

scatterplot(xup)
grid on
title('Constelación de simbolos generados');
set(gcf, 'Position', [50 50 600 600],'Color', 'w');
