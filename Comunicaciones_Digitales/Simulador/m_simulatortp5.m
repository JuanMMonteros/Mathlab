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
    config_s.tx_s.M = 4;                       % Cantidad de niveles de la modulacion
    config_s.tx_s.NOS = 4;                      % Tasa de sobremuestreo
    config_s.tx_s.Lsymbs = 1e6;                 % Cantidad de simbolos
    config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
    config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
    config_s.tx_s.pulse_shaping_type = 0;       % 0: RRC, 1: RC
    %ch
    config_s.ch_awgn.EbNo_db =3; 
    config_s.ch_awgn.ISE = 1; %1 activada, 0 desacticada
    config_s.ch_awgn.firorder = 17;
    config_s.ch_awgn.fc = 16e9;
    %AGC
    config_s.agc.target = 0.3;
    config_s.agc.adc_phase = 1;
    %Ecualizador Adaptivo
    config_s.ec_s.ntap = 63; 
    config_s.ec_s.N =2;        
    config_s.ec_s.step_cma=2^-11;
    config_s.ec_s.step_dd=2^-11;
    config_s.ec_s.tap_leak_gain=1e-10;
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
    config_s.ber_s.phase=config_s.agc.adc_phase;
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
    FSE=o_ec_s.W;
    Chanel_F=ch_awgn.b;
    H=o_tx_s.filter;
    
%Respuesta al impulso del FFE
figure %1.A
subplot 211
stem(real(FSE))
title('Parte Real')
subplot 212
stem(imag(FSE))
title('Parte Imaginaria')
 %1.b Respues en Frecuencia del canal y del ecualizador y la convolucion de los
 %mismo 
     figure %b

      NFFT=1024;
      fs=config_s.tx_s.BR*config_s.tx_s.NOS;
      CH =fftshift(fft(Chanel_F, NFFT));
      CH=CH/CH(NFFT/2);
      f=[-fs/2:fs/NFFT:fs/2-fs/(NFFT)];
      plot(f/1e9, 20*log10(abs(CH)), '--k', 'LineWidth',2);
      grid on
      hold on
      EC =fftshift(fft(FSE, NFFT));
      EC=EC/EC(NFFT/2);
      plot(f/1e9, 20*log10(abs(EC)),'-r', 'Linewidth',2);
      A = fftshift(fft(conv(FSE,Chanel_F), NFFT));
      A=A/A(NFFT/2);
      plot(f/1e9, 20*log10(abs(A)),'-g', 'Linewidth',2);
       
       xlabel('Frequency [GHz]')
       ylabel('PSD Magnitude ')
       ylim([-100,20]);
       title('Respuesta En Frecuencia')
       set(gcf, 'Position', [50 50 500 500],'Color', 'w');
       legend({'Respuesta del canal','Respusta del Ecualizador','Respuesta de La convolucion'},'Location','s') 
    end
    
    o_data_s.ber_theo = ber_check.ber_theo;
    o_data_s.ber_sim = ber_check.ber_sim ;

end