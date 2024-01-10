%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Test Unitario para el receptor                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
close all; 
clear

%% -- SETTINGS -- %

rx_config.BR = 32e9;       % Baud rate
rx_config.M = 16;          % Cantidad de niveles de la modulacion
rx_config.NOS = 2;         % Tasa de sobremuestreo
rx_config.dwn_phase = 0;   % fase de downsampling

% -- settings TX -- %
M = rx_config.M;    % Order modul
BR = rx_config.BR;  % Baud Rate
N = rx_config.NOS;  % Oversampling rate
EbNo_db = 8;        % EbNo in dB
L = 1e5;            % Simulation Length
rolloff = 0.5;      % Pulse shaping rolloff
h_taps = 201;       % Pulse shaping taps

% -- settings Ch -- %   
bw  = 20e9;         % Bandwid
bw_order = 18;      % Bw order    
ebno_db = EbNo_db;  % EbNo in dB
enable_noise = 0;   % noise
en_bw_lim = 1;      % Bw limit


                   % -- basic Tx -- %

fs = N*BR;  % Sampling rate to emulate analog domain
T = 1/BR;   % Time interval between two consecutive symbols
Ts = 1/fs;  % Time between 2 consecutive samples at Tx output

% generate labels
dec_labels = randi([0 M-1], L, 1); 
% generate simbol
ak = qammod(dec_labels,M); 
% Upsampling to change sampling rate
xup = N * upsample(ak,N);
% Filter to interpolate the signal
htx_v = RootRainCosine(BR/2, fs, rolloff, h_taps, 0);
% output filter tx
yup = filter(htx_v,1,xup);


                   % -- basic Ch -- %
if en_bw_lim == 1 
    fs = NOS*BR;             % frec sampling
    fc = bw/(fs/2);          % Frec de corte
    hch_v = fir1(bw_order,fc); % filtro de canal     
else
    hch_v = 1;
end      

% -- Generación de EbNo -- %
k = log2(M);    
ebno_veces = 10^(ebno_db/10);
SNR_slc = k*ebno_veces;
SNR_ch = SNR_slc/NOS;
Pot_sign = var(yup);
Pot_noise = Pot_sign / SNR_ch;  

% Ruido agregado de ruido al canal 
if enable_noise == 1
    noise_v = sqrt(Pot_noise/2).*(randn(size(yup)) + 1j.*randn(size(yup)));
else
    noise_v = 0;
end
% -- Filtrado del canal -- %
ych_v = yup + noise_v;  % Add noise
ych_v = filter(hch_v, 1, ych_v); % Out filter ch

% -- PROCESS -- %
i_rx_s.htx_v = htx_v; % Filtro tx    
i_rx_s.hch_v = hch_v; % Filtro ch 
i_rx_s.signal_v = ych_v; % Signal out ch

o_rx_s = ReceptorMQAM(i_rx_s, rx_config);

% -- Plots -- %

scatterplot(o_rx_s.y_v);
title('Constelación entrada del slicer');
set(gcf, 'Position', [50 50 600 600],'Color', 'w');
