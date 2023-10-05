function [o_data_s] = BER_CHECKER(i_config_s,rx,tx_symbs)
%--------------------------%
%     DEFAULT SETTINGS
 %--------------------------%
 
    config.M = 4;  % Cantidad de niveles de la modulacion
    config.EbNo_db = 10;
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
    M = config.M;
    EbNo_db=config.EbNo_db;
    %--------------------------%
    %         PROCESS
    %--------------------------%
    % Slicer
ak_hat = my_slicer(rx, M);
    % salidas del slicer

% Theo ber
ber_theo = berawgn(EbNo_db, 'qam', M);

% Estimated ber
% se recomiendo simular hasta contar al menos 100 errores
[ber_sim, n_errors] = my_ber_checker(ak_hat, tx_symbs/sqrt(1), M, 'auto');

o_data_s.ber_theo=ber_theo;
o_data_s.ber_sim=ber_sim;

end

