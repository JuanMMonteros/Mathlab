function [o_data_s] = Chanel_AWGN_IBN(i_config_s,s)
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
    
    f=(BR+rolloff*BR)/2*NOS;
    fc=fc/f;
    if ISE
        b=fir1(firorder,fc);
    else
        b=1;
    end
    s = filter(b,1,s);
    
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
o_data_s.yup_n=rx; 
o_data_s.b=b;

end

