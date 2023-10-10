%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Juan Monteros And Matias Mollcker
% Created on   : 2023
% Description  : QAM simulator
%-----------------------------------------------------------------------------%
 
function o_data_s = m_simulatortp5(i_cfg_s)
close all 
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    % -- General --
    config_s.en_plots = 1;

    % -- Tx --
    config_s.tx_s.BR = 32e9;                    % Baud rate
    config_s.tx_s.M = 16;                       % Cantidad de niveles de la modulacion
    config_s.tx_s.NOS = 4;                      % Tasa de sobremuestreo
    config_s.tx_s.Lsymbs = 1e6;                 % Cantidad de simbolos
    config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
    config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
    config_s.tx_s.pulse_shaping_type = 0;       % 0: RRC, 1: RC
    %ch
    config_s.ch_awgn.EbNo_db = 14; 
    config_s.ch_awgn.ISE = 1; %1 activada, 0 desacticada
    config_s.ch_awgn.firorder = 17;
    config_s.ch_awgn.fc = 20e9;
    %AGC
    config_s.agc.target = 0.3;
    config_s.agc.adc_phase = 1;
    %Ecualizador Adaptivo
    config_s.ec_s.ntap = 63; 
    config_s.ec_s.N =2;        
    config_s.ec_s.step_cma=2^-11;
    config_s.ec_s.step_dd=2^-11;
    config_s.ec_s.tap_leak_gain=1e-4;
    config_s.ec_s.force_cma_enable=0;
    
    

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
    config_s.rx_s.ntaps = config_s.tx_s.pulse_shaping_ntaps;
    config_s.ber_s.M = config_s.tx_s.M;
    config_s.ber_s.EbNo_db = config_s.ch_awgn.EbNo_db;
    config_s.agc.NOS =  config_s.tx_s.NOS;
    config_s.agc.N =  config_s.ec_s.N;
    config_s.ec_s.M =  config_s.tx_s.M; 

    %--------------------------%
    %          PROCESS
    %--------------------------%

    % -- Tx --
    o_tx_s = transmisor_MQAM(config_s.tx_s);
    % -- CH --
    %ch_awgn =  Chanel_AWGN_IBN(config_s.ch_awgn,o_tx_s.oversampled_output);
    ch_awgn =  Chanel_AWGN(config_s.ch_awgn,o_tx_s.oversampled_output);
    %-- AGC --
    agc = AGC(config_s.agc,ch_awgn.yup_n);
    % -- Rx --
    o_ec_s = Ecualizador(config_s.ec_s,agc.rx_norm,o_tx_s.RCMA);
    % -- slicer y berchequer
    
    ber_check = BER_CHECKER(config_s.ber_s,o_ec_s.yk,o_tx_s.tx_symbs);
    
    if config_s.en_plots 

       tx_data = o_tx_s.oversampled_output;
       y=ch_awgn.yup_n;
       r=o_ec_s.yk;

       
       eyediagram(tx_data(13:23e3), config_s.tx_s.NOS );%a
       
       figure %b
       NFFT = 1024;
       WELCH_OVERLAP = 0;
       fs=config_s.tx_s.BR*config_s.tx_s.NOS;
       [Pxx, f] = pwelch(tx_data, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs);
       Px_dB= 10*log10(Pxx);
       Px_dB = Px_dB - Px_dB(1); %Normalize Px_dB to get 0dB at f=0, only for plot purposes
       plot(f/1e9, Px_dB,'-b', 'Linewidth',2)
       grid on
       hold on 
       [Pxx, f] = pwelch(y, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs);
       Px_dB= 10*log10(Pxx);
       Px_dB = Px_dB - Px_dB(1); %Normalize Px_dB to get 0dB at f=0, only for plot purposes
       plot(f/1e9, Px_dB,'-r', 'Linewidth',2)
       xlabel('Frequency [GHz]')
       ylabel('PSD Magnitude ')
       title('PSD salida Tx vs Entrada RX')
       set(gcf, 'Position', [50 50 500 500],'Color', 'w');
       legend({'Salida del Transmisor','Entrada del receptor'},'Location','s')

       figure %c
       [Pxx, f] = pwelch(y, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs);
       Px_dB= 10*log10(Pxx);
       Px_dB = Px_dB - Px_dB(1); %Normalize Px_dB to get 0dB at f=0, only for plot purposes
       plot(f/1e9, Px_dB,'-b', 'Linewidth',2)
       grid on
       hold on 
       [Pxx, f] = pwelch(r, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs);
       Px_dB= 10*log10(Pxx);
       Px_dB = Px_dB - Px_dB(1); %Normalize Px_dB to get 0dB at f=0, only for plot purposes
       plot(f/1e9, Px_dB,'-r', 'Linewidth',2)
       xlabel('Frequency [GHz]')
       ylabel('PSD Magnitude ')
       title('PSD Entrada Rx vs Salida RX')
       set(gcf, 'Position', [50 50 500 500],'Color', 'w');
       legend({'Entrada del Reseptor','Salida del Receptor'},'Location','s')
        
      scatterplot(r(100e3:101e3)); %d
      
      figure; %e
      histogram(real(r), 'BinMethod', 'auto');
      title('Histograma de Parte Real antes del Slicer');
      xlabel('Valores');
      ylabel('Frecuencia');
      figure;
      histogram(imag(r), 'BinMethod', 'auto');
      title('Histograma de Parte Imaginaria antes del Slicer');
      xlabel('Valores');
      ylabel('Frecuencia');
        

    end
    
    o_data_s.ber_theo = ber_check.ber_theo;
    o_data_s.ber_sim = ber_check.ber_sim ;

end