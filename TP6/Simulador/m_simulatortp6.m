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
    position_and_size = [100 100 800 600];  
    BR=32e9;
    

    % -- General --
    config_s.en_plots = 0;

    % -- Tx --
    config_s.tx_s.BR = 32e9;                    % Baud rate
    config_s.tx_s.M = 16;                       % Cantidad de niveles de la modulacion
    config_s.tx_s.NOS = 1;                      % Tasa de sobremuestreo
    config_s.tx_s.Lsymbs = 1e6;                 % Cantidad de simbolos
    config_s.tx_s.rolloff = 0.5;                % Rolloff del filtro conformador
    config_s.tx_s.pulse_shaping_ntaps = 201;    % Cantidad de taps del PS
    config_s.tx_s.pulse_shaping_type = 0;       % 0: RRC, 1: RC
    %---ch---
    config_s.ch_awgn.EbNo_db =10; 
    config_s.ch_awgn.ISE = 0; %1 activada, 0 desacticada
    config_s.ch_awgn.firorder = 17;
    config_s.ch_awgn.fc = 20e9;
    %ch Configracion de portadora 
    config_s.ch_awgn.delta_freq = 10e6; % Offset del LO
    config_s.ch_awgn.LW = 00e3; % Ancho de linea [Hz]
    config_s.ch_awgn.theta0 = 15/180*pi; % Fase inicial
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
    ber_check = BER_CHECKER(config_s.ber_s,o_tx_s.tx_symbs,rx_s.pll_output,rx_s.ak_hat);
    % -- Ch -- %
    ch_out_symb = ch_awgn.yup_n;
    phase_noise =  ch_awgn.phase_noise;
    
    
    % -- Rx Carrie Recovery -- %
    pll_output = rx_s.pll_output;      
    nco_output = rx_s.nco_output;
    integral_branch = rx_s.integral_branch; 
    ak_hat = rx_s.ak_hat;
    phase_error = rx_s.phase_error;
    Ldata =  rx_s.Ldata;
    fz=16;

    %Constelacion Entrada PLL
    scatterplot(ch_out_symb(1000:10000))
    grid on 
    legend("Entrad Escalon  de 45 grados", 'Location', 'eastoutside', 'FontSize', fz-5);
    title("Constelacion a la Entrada del PLL");
    set(gcf, 'Position', position_and_size, 'Color', 'w');
    %Constelacion Salida PLL
    scatterplot(pll_output(Ldata/2:(Ldata/2+9000)))
    grid on;
    title("Contelacion a la Salida del PLL");
    
    
    % Calcular la densidad espectral de potencia (PSD) usando pwelch
    fs = BR; % Frecuencia de muestreo
    window = hamming(length(rx_s.nco_output)); % Ventana (puedes cambiar a otra si prefieres)
    noverlap = 0; % Sin superposici�n (puedes ajustar esto seg�n sea necesario)
    nfft = 1024; % Por defecto utiliza el pr�ximo valor de potencia de 2 de la longitud de la se�al

    [Pxx, f] = pwelch(rx_s.nco_output, window, noverlap, nfft, fs);
    % Obtener la potencia del tono recuperado (primera frecuencia)
    Px_dB= 10*log10(Pxx);
    Px_dB = Px_dB - Px_dB(1);
   % Graficar la densidad espectral de potencia (PSD)
       grid on
    figure;
   plot(f, Px_dB,'-b', 'Linewidth',2)
   %xlim([10e7:10e10]);
   title('Densidad Espectral de Potencia (PSD)');
   grid on
   xlabel('Frecuencia (Hz)');
   ylabel('PSD (dB/Hz)');
        if config_s.en_plots
        % -- Ruido de fase del Ch -- %
        figure
        plot(phase_noise(1000:10000))
        grid on
        tit = sprintf('Ruido de fase del canal');
        name = sprintf('PhaseNoise');
        legend(name,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
        title(tit,'Interpreter','latex', 'FontSize', fz);
        xlabel('Time')
        ylabel('Channel Phase Noise [rad]')
        set(gcf, 'Position',position_and_size,'Color', 'w');
        
%%
        % -- Error de fase del receptor -- %
        figure
        plot(phase_error(1000:10000))
        grid on
        tit = sprintf('Error de fase');
        name = sprintf('PhaseError');
        legend(name,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
        title(tit,'Interpreter','latex', 'FontSize', fz);        
        xlabel("Time")
        ylabel("Phase Error")
        set(gcf, 'Position',position_and_size,'Color', 'w');
        
    
        % -- Rama integral del PLL -- %
        figure
        plot(integral_branch(1000:10000)*BR/2/pi/1e6)
        grid on
        tit = sprintf('Rama integral del PLL');
        name = sprintf('IntegralBranch');
        legend(name,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
        title(tit,'Interpreter','latex', 'FontSize', fz);          
        xlabel("Time")
        ylabel("Integral Branch [MHz]")
        set(gcf, 'Position',position_and_size,'Color', 'w');
        

        % -- Salida del NCO vs error de fase -- %
        figure
        plot(nco_output(1000:10000)); hold all
        plot(phase_noise(1000:10000)); hold all
        grid on
        tit = sprintf('Salida del NCO vs error de fase');
        name1 = sprintf('PLL_NCO');
        name2 = sprintf('Channel');
        legend(name1,name2,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
        title(tit,'Interpreter','latex', 'FontSize', fz);        
        xlabel("Time")
        ylabel("Phase [rad]")
        set(gcf, 'Position',position_and_size,'Color', 'w');
       
        
        scatterplot(pll_output(Ldata/2:(Ldata/2+9000)))
        grid on;
        tit = sprintf('Salida del PLL');
        if config_s.tx_s.M == 4
            name = sprintf('QAM');
        elseif config_s.tx_s.M > 4
            name = sprintf('QAM%d',config_s.tx_s.M);
        end
        legend(name,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
        title(tit,'Interpreter','latex', 'FontSize', fz);
        set(gcf, 'Position',position_and_size,'Color', 'w');
       
        
        % -- Entrada y salida del slicer -- %
        figure
        plot(real(pll_output(1000:10000)), '.r');
        hold all
        plot(real(ak_hat(1000:10000)), 'xb');
        grid on;
        tit = sprintf('Entrada y salida del slicer');        
        name1 = sprintf('Slicer Input');
        name2 = sprintf('Slicer Output');
        legend(name1,name2,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
        title(tit,'Interpreter','latex', 'FontSize', fz);
        set(gcf, 'Position',position_and_size,'Color', 'w');      
    end
    
    o_data_s.ber_theo = ber_check.ber_theo;
    o_data_s.ber_sim_sin_corregir = ber_check.ber_sim_sin_corregir;
    o_data_s.ber_sim_corrigiendo = ber_check.ber_sim_corrigiendo;

end
