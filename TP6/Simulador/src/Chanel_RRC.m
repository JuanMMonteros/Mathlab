function [o_data_s] = Chanel_RRC(i_config_s,s)
%CHANEL_AWGN Summary of this function goes here
%--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    config.EbNo_db = 5; 
    config.M = 4;   % Cantidad de niveles de la modulacion
    config.NOS = 2;
    config.ISE=0;
    config.firorder=1;
    config.fc=20e9;
    config.BR=32e9;
    config.rolloff=0.5;
    config.delta_freq = 0e6; % Offset del LO
    config.LW = 00e3; % Ancho de linea [Hz]
    config.theta0 = 0/180*pi; % Fase inicial
    config.frequency_fluctuations_amp = 0e6;
    config.frequency_fluctuations_freq = 1e3;
    config.phase_tone_amp = 10/180*pi;
    config.phase_tone_freq = 250e6;

    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%

    fn = fieldnames(i_config_s);
    for k = 1:numel(fn)
        if isfield(config,(fn{k}))==1
            config.(fn{k})= i_config_s.(fn{k});
        else
            error("%s: Parametro del simulador no valido", fn{k})
        end
    end

    %--------------------------%
    %         VARIABLES
    %--------------------------%
    
    EbNo_db = config.EbNo_db;
    M = config.M;
    NOS = config.NOS;
    ISE = config.ISE;
    BR=config.BR;
    fc= config.fc;
    firorder= config.firorder;
    rolloff=config.rolloff;
    delta_freq =  config.LW; % Offset del LO
    LW = 00e3; % Ancho de linea [Hz]
    theta0 = config.theta0; % Fase inicial
    frequency_fluctuations_amp = config.frequency_fluctuations_amp;
    frequency_fluctuations_freq = config.frequency_fluctuations_freq;
    phase_tone_amp = config.phase_tone_amp;
    phase_tone_freq = config.phase_tone_freq;

    
    f=(BR+rolloff*BR)/2*NOS;
    fc=fc/f;
    Ts=1/BR;
    fs=1/Ts;
    
    %--------------------------%
    %         PROCESS
    %--------------------------%


% EbNo to channel snr
k = log2(M);
EbNo = 10^(EbNo_db/10);
SNR_slc = EbNo * k;
SNR_ch = SNR_slc / NOS;

% Noise power
Ps = var(s);
Pn = Ps/SNR_ch;

% Noise generator
n = sqrt(Pn/2) .* (randn(length(s),1) + 1j.*randn(length(s),1));

% Noise addition
rx = s + n; 
 if ISE
         b=fir1(firorder,fc);
         rx_noisy = filter(b,1,rx);
     else
         rx_noisy = rx;
 end
     

% Agregar efectos de portadora
Ldata= length(rx_noisy);
time = (0:Ldata-1).'.*Ts;

lo_offset = exp(1j*2*pi*delta_freq*time); % LO OFFSET
phase_offset = exp(1j*theta0); % STATIC PHASE
% Fluctuaciones de frecuencia
freq_fluctuations = exp(1j*frequency_fluctuations_amp/frequency_fluctuations_freq.*sin(2*pi*frequency_fluctuations_freq.*time));
% Para el ruidod de fase
freq_noise = sqrt(2*pi*LW/fs).*randn(Ldata,1);
phase_noise = cumsum(freq_noise); % Proceso de Wiener
osc_pn = exp(1j.*phase_noise);
phase_tone = exp(1j.*phase_tone_amp.*sin(2*pi*phase_tone_freq.*time));
rxs = rx_noisy.*lo_offset.*phase_offset.*osc_pn.*freq_fluctuations.*phase_tone;
%rxs =rx_noisy.*lo_offset.*phase_offset.*osc_pn;

o_data_s.yup_n=rxs;
o_data_s.phase_noise = phase_noise;

end