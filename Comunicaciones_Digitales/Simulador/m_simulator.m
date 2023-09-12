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
    config_s.tx_s.M = 4;                       % Cantidad de niveles de la modulacion
    config_s.tx_s.NOS = 2;                      % Tasa de sobremuestreo
    config_s.tx_s.Lsymbs = 1e5;                 % Cantidad de simbolos
    config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
    config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
    config_s.tx_s.pulse_shaping_type = 1;       % 0: RRC, 1: RC
    config_s.ch_awgn.EdNo_db = 10; 
     config_s.ch_awgn.M =  config_s.tx_s.M;                       % Cantidad de niveles de la modulacion
    config_s.ch_awgn.NOS =  config_s.tx_s.NOS;
    config_s.rx_s.filter_type = 1 ;
    config_s.rx_s.NOS =  config_s.tx_s.NOS;
    config_s.rx_s.ntaps = config_s.tx_s.pulse_shaping_ntaps; 
       
    

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
    % -- CH --
    ch_awgn =  Chanel_AWGN(config_s.ch_awgn,o_tx_s.oversampled_output);
    % -- Rx --
    o_rx_s = Reseptor(config_s.rx_s,ch_awgn.yup_n,o_tx_s.filter);

    %--------------------------%
    %          PLOTS
    %--------------------------%

    if config_s.en_plots 

        rx_symbs = o_rx_s.output_rx;
        tx_data = ch_awgn.yup_n;

        scatterplot(rx_symbs);
        histogram(real(tx_data), 'BinMethod', 'auto');
        eyediagram(tx_data(1:10e3), config_s.tx_s.NOS );

        figure
        pwelch(tx_data)

    end

    %--------------------------%
    %          OUTPUT
    %--------------------------%

    o_data_s.tx_sym_v = o_tx_s.tx_symbs;
    o_data_s.ber_theo = config_s.tx_s.M ;

end