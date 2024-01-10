%-------------------------------------------------------------------------%
%                                 FULGOR
% Programmer(s): Daniel A. Juarez
% Created on   : October 2023
% Description  : Transmission channel Ch
%-------------------------------------------------------------------------%

function [o_data_s] = CanalDeTransmision(i_ch_s, i_config_s)

    %----------------------------------% 
    %         DEFUALT SETTINGS
    %----------------------------------%
                % -- Parametros de Portadora -- %
    config.delta_freq = 0e6;                    % Offset del LO
    config.LW = 00e3;                           % Ancho de linea [Hz]
    config.theta0 = 0/180*pi;                   % Fase inicial
    config.frequency_fluctuations_amp = 0e6;    % Fluc de amplitud
    config.frequency_fluctuations_freq = 1e3;   % Fluct de Frec
    config.phase_tone_amp = 10/180*pi;          % Amplitud del tono
    config.phase_tone_freq = 250e6;             % Fase del tono
                % -- Parametros del canal -- %
    config.EbNo_dB = 6;                         % EbNo en [dB]
    config.fs = 32e3;                           % Frecuencia de sampling
    config.M  = 4;                              % Nivel de modulación  
    config.Ts = 1/32e3;                         % Tiempo de sampling 

    %----------------------------------% 
    %           REASSIGNMENT
    %----------------------------------%
    % Sobre escritura de los parámetros
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
    rand('seed',1);             % Semillas rand
    randn('seed',51);

    % -- Parámetros de entrada -- %
    xsymb = i_ch_s.xsymb;       
    
    % -- Parámetros de portadora -- %
    delta_freq  = config.delta_freq;                                   
    LW  = config.LW;                                                   
    theta0  = config.theta0;                                           
    frequency_fluctuations_amp  = config.frequency_fluctuations_amp;   
    frequency_fluctuations_freq = config.frequency_fluctuations_freq;  
    phase_tone_amp  = config.phase_tone_amp;                           
    phase_tone_freq = config.phase_tone_freq;   

    % -- Parámetros del canal -- %
    M = config.M;              
    EbNo_dB = config.EbNo_dB;   
    Ts = config.Ts;             
    fs = config.fs;

    %--------------------------%
    %         PROCESS
    %--------------------------%
    
    % -- Calculos de la SNR -- %
    EbNo = 10^(EbNo_dB/10);
    SNR_Slc = EbNo * log2(M);
    SNR_canal = SNR_Slc; 
    
    % -- Agregado de ruido al canal -- %
    Psignal = var(xsymb);
    Pnoise = Psignal/SNR_canal;
    No = Pnoise;                                % varianza por unidad de frecuencia
    noise_real = sqrt(No/2)*randn(size(xsymb)); % Lo que acompaña al ruido es el sigma (desvio estandard)
    noise_imag = sqrt(No/2)*randn(size(xsymb));
    noise = noise_real + 1j*noise_imag;
    % NO VALE noise= (1+1j)*randn(size(yup))
    rx_noisy = xsymb + noise;
    clear noise noise_real noise_imag

    % -- Rotacion de portadora (Agregado de efectos de portadora) -- % 
    Ldata= length(rx_noisy);
    time = (0:Ldata-1).'.*Ts;

    lo_offset = exp(1j*2*pi*delta_freq*time); % LO Offset
    phase_offset = exp(1j*theta0); % Static Phase

    %  -- Fluctuaciones de frecuencia -- %
    % freq_fluctuations = exp(1j*frequency_fluctuations_amp/ ...
    % frequency_fluctuations_freq.*sin(2*pi*frequency_fluctuations_freq.*time));

    % -- Para el ruido de fase -- %
    freq_noise = sqrt(2*pi*LW/fs).*randn(Ldata,1);
    phase_noise = cumsum(freq_noise); % Proceso de Wiener
    osc_pn = exp(1j.*phase_noise);
    % phase_tone = exp(1j.*phase_tone_amp.*sin(2*pi*phase_tone_freq.*time));
    % rxs = rx_noisy.*lo_offset.*phase_offset.*osc_pn.*freq_fluctuations.*phase_tone;
    ch_out_symb = rx_noisy.*lo_offset.*phase_offset.*osc_pn;

    %--------------------------%
    %         OUTPUT
    %--------------------------%
    
    o_data_s.ch_out_symb = ch_out_symb; % Signal out filter ch
    o_data_s.phase_noise = phase_noise; % Ruido de fase
end    
