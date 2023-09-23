
function [o_data_s] = Reseptor(i_config_s,r,h)
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    config.filter_type = 1; 
    config.NOS = 2;  % Cantidad de niveles de la modulacion
    config.ntaps=201;
    
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
    
    filter_type= config.filter_type;
    NOS = config.NOS;
    n_taps=config.ntaps;
    %--------------------------%
    %         PROCESS
    %--------------------------%
%% RX
% Filter with MF
switch filter_type
    case 1      
         f = conj(h(end:-1:1)); % Filter with MF
    case 2
        f=1;
end
h_delay = (n_taps-1)/2;
ymf = filter(f,1,[r; zeros(h_delay, 1)]);
ymf = ymf(1+h_delay:end);

% Downsampling
rx_down = ymf(1:NOS:end);

o_data_s.output_rx= rx_down;

end

