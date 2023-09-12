function [o_data_s] = Chanel_AWGN(i_config_s,s)
%CHANEL_AWGN Summary of this function goes here
%--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    config.EdNo_db = 5; 
    config.M = 4;   % Cantidad de niveles de la modulacion
    config.NOS = 2;  % Cantidad de niveles de la modulacion
    
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
    
    EdNo_db = config.EdNo_db;
    M = config.M;
    NOS = config.NOS;
    %--------------------------%
    %         PROCESS
    %--------------------------%


% EbNo to channel snr
k = log2(M);
EbNo = 10^(EdNo_db/10);
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

end

