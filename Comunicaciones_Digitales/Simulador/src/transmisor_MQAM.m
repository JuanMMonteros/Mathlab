function [o_data_s] = transmisor_MQAM(i_config_s)
    
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    config.BR = 32e9;                   % Cantidad de niveles de la modulacion
    config.M = 4;                       % Cantidad de niveles de la modulacion
    config.NOS = 2;                     % Tasa de sobremuestreo
    config.Lsymbs = 100e3;              % Cantidad de simbolos
    config.rolloff = 0.1;               % Rolloff del filtro conformador
    config.pulse_shaping_ntaps = 201;   % Taps del filtro conformador
    config.pulse_shaping_type = 0;      % 0: RRC, 1: RC
  
    
   
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
    
    BR = config.BR;
    M = config.M;
    NOS = config.NOS;
    Lsymbs = config.Lsymbs;
    rolloff = config.rolloff;
    pulse_shaping_ntaps = config.pulse_shaping_ntaps;
    pulse_shaping_type = config.pulse_shaping_type;
    
    fs = BR*NOS;

    %--------------------------%
    %         PROCESS
    %--------------------------%

    % QAM Symbols generation 
    dec_labels = randi([0 M-1], Lsymbs, 1);
    tx_symbs = sqrt(1)*qammod(dec_labels,M);

    % Upsampling to change sampling rate
    xup = upsample(tx_symbs, NOS);

    % Pulse shaping con RRC-RC
    if pulse_shaping_type==0
        htx = root_raised_cosine(BR/2, fs, rolloff, pulse_shaping_ntaps, 0);
    else
        htx = raised_cosine(BR/2, fs, rolloff, pulse_shaping_ntaps, 0);
    end
 
    %s = filter(htx,1,xup);
h_delay = (pulse_shaping_ntaps-1)/2;
s = filter(htx,1,[xup; zeros(h_delay, 1)]);
s =NOS*s(1+h_delay:end);

    clear xup
    
    %--------------------------%
    %         OUTPUT
    %--------------------------%
    RCMA = sqrt(mean(abs(tx_symbs(1:1000)).^4)/mean(abs(tx_symbs(1:1000)).^2));

    o_data_s.oversampled_output = s;
    o_data_s.filter = htx;
    o_data_s.tx_symbs = tx_symbs;
    o_data_s.RCMA = RCMA;
end

