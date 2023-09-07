%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero 
% Created on   : 2023
% Description  : QAM simulator
%-----------------------------------------------------------------------------%

function o_data_s = m_simulator(i_cfg_s)

    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    % -- General --
    config_s.en_plots = 1;

    % -- Tx --
    config_s.tx_s.BR = 32e9;                    % Baud rate
    config_s.tx_s.M = 16;                       % Cantidad de niveles de la modulacion
    config_s.tx_s.NOS = 2;                      % Tasa de sobremuestreo
    config_s.tx_s.Lsymbs = 1e5;                 % Cantidad de simbolos
    config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
    config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
    config_s.tx_s.pulse_shaping_type = 1;       % 0: RRC, 1: RC
       
    

    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%

    if nargin > 0
        config_s = overwrite_parameters(i_cfg_s, config_s);
    end

    %--------------------------%
    %          PROCESS
    %--------------------------%

    % -- Tx --
    o_tx_s = transmisor_MQAM(config_s.tx_s);

    % -- Rx --

    %--------------------------%
    %          PLOTS
    %--------------------------%

    if config_s.en_plots 

        tx_symbs = o_tx_s.tx_symbs;
        tx_data = o_tx_s.oversampled_output;

        scatterplot(tx_symbs)
        eyediagram(tx_data(1:10e3), config_s.tx_s.NOS )

        figure
        pwelch(tx_data)

    end

    %--------------------------%
    %          OUTPUT
    %--------------------------%

    o_data_s.tx_sym_v = o_tx_s.tx_symbs;
    o_data_s.ber_theo = config_s.tx_s.M ;

end
