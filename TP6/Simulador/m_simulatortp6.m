%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Juan Monteros And Matias Mollecker
% Created on   : 2024
% Description  : RRC simulator
%-----------------------------------------------------------------------------%
 
function o_data_s = m_simulatortp6(i_cfg_s)
close all 
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%                                 
    

    % -- General --
    config_s.en_plots = 0;

    % -- Tx --
    config_s.tx_s.BR = 32e9;                    % Baud rate
    config_s.tx_s.M = 16;                       % Cantidad de niveles de la modulacion
    config_s.tx_s.NOS = 1;                      % Tasa de sobremuestreo
    config_s.tx_s.Lsymbs = 1e4;                 % Cantidad de simbolos
    config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
    config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
    config_s.tx_s.pulse_shaping_type = 0;       % 0: RRC, 1: RC
    %---ch---
    config_s.ch_awgn.EbNo_db =10; 
    config_s.ch_awgn.ISE = 0; %1 activada, 0 desacticada
    config_s.ch_awgn.firorder = 17;
    config_s.ch_awgn.fc = 20e9;
    %ch Configracion de portadora 
    config_s.ch_awgn.delta_freq = 50e6; % Offset del LO
    config_s.ch_awgn.LW = 00e3; % Ancho de linea [Hz]
    config_s.ch_awgn.theta0 = 0/180*pi; % Fase inicial
    config_s.ch_awgn.frequency_fluctuations_amp = 0e6;
    config_s.ch_awgn.frequency_fluctuations_freq = 1e3;
    config_s.ch_awgn.phase_tone_freq = 250e6;
    %--PLL--
    config_s.rx_s.Kp =0.05;
    config_s.rx_s.Ki = config_s.rx_s.Kp/1000;
    

    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%

    if nargin > 0
        config_s = overwrite_parameters(i_cfg_s, config_s);
    end
    %Sharrimg paramiter
    config_s.ch_awgn.M =  config_s.tx_s.M;                       % Cantidad de niveles de la modulacion
    config_s.ch_awgn.NOS =  config_s.tx_s.NOS;
    config_s.ch_awgn.BR =config_s.tx_s.BR;
    config_s.ch_awgn.rolloff=config_s.tx_s.rolloff;
    config_s.rx_s.NOS =  config_s.tx_s.NOS;
    config_s.rx_s.M =  config_s.tx_s.M;
    config_s.ber_s.M = config_s.tx_s.M;
    config_s.ber_s.EbNo_db = config_s.ch_awgn.EbNo_db;
    config_s.ber_.Ldata = config_s.tx_s.Lsymbs;

    %--------------------------%
    %          PROCESS
    %--------------------------%

    % -- Tx --
    o_tx_s = transmisor_MQAM(config_s.tx_s);
    % -- CH --
    %ch_awgn =  Chanel_AWGN_IBN(config_s.ch_awgn,o_tx_s.oversampled_output);
    ch_awgn =  Chanel_RRC(config_s.ch_awgn,o_tx_s.tx_symbs);
    %--PLL RX--
    rx_s = RX_RCC(config_s.rx_s,ch_awgn.yup_n);
    %--Ber Cheker--
    % ber_check = BER_CHECKER(config_s.ber_s,o_tx_s.tx_symbs,rx_s.pll_output,rx_s.ak_hat);
    % -- Ch -- %
    ch_out_symb = ch_awgn.yup_n;    
    % -- Rx Carrie Recovery -- %
    pll_output = rx_s.pll_output;
    Ldata =  rx_s.Ldata;

    if config_s.en_plots
    %Constelacion Entrada PLL
    scatterplot(ch_out_symb(1000:10000))
    grid on 
    title("Constelacion a la Entrada del PLL");
    %Constelacion Salida PLL
    scatterplot(pll_output(Ldata/2:(Ldata/2+4000)))
    grid on;
    title("Contelacion a la Salida del PLL");
    end
        
    
    %o_data_s.ber_theo = ber_check.ber_theo;
    %o_data_s.ber_sim_sin_corregir = ber_check.ber_sim_sin_corregir;
    %o_data_s.ber_sim_corrigiendo = ber_check.ber_sim_corrigiendo;
    o_data_s.pll_output=rx_s.pll_output;
    o_data_s.tita_in=ch_awgn.phase_tone;

end
