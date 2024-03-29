%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero 
% Created on   : September 2023
% Description  : BPS simulator
%-----------------------------------------------------------------------------%

function o_data_s = bps_sim(i_cfg_s)

    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%
    %-- plots --$
    cfg_s.plots=1;
    % -- Tx --
    cfg_s.tx_s.n_sym                  = 1e5     ;   % Number of symbols
    cfg_s.tx_s.n_pil                  = 0*100     ;   % Number of pilots symbols
    cfg_s.tx_s.n_payload              = 127     ;   % Symbols between pilots
    
    cfg_s.tx_s.M                      = 16      ;   % Modulation levels
    cfg_s.tx_s.BR                     = 32e9   ;   % Baud rate [Bd]
    
    % -- Laser --
    cfg_s.lw                          = 10e3       ;   % LW [Hz]
    
    % -- Noise --
    cfg_s.ebno_db                     = 6.7      ;   % EbNo [dB] 
    
    % -- BPS --
    cfg_s.bps_s.n_phases              = 4     ;   % BPS phases
    cfg_s.bps_s.filter_length         = 101   ;   % BPS filter len
    
    % -- CS corrector --
    cfg_s.en_cs_corrector             = 1       ;   % CS corrector for BPS path
    cfg_s.cs_cor_s.window_length      = 64      ;   % CS corrector window len
    
    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%

    if nargin > 0
        cfg_s = overwrite_parameters(i_cfg_s, cfg_s);
    end

    clear i_cfg_s;
    
    % BPS
    cfg_s.bps_s.mapper_levels = cfg_s.tx_s.M;
    
    % Round filters len to odd
    cfg_s.bps_s.filter_length = round_odd(cfg_s.bps_s.filter_length);
    
    %--------------------------%
    %          PROCESS
    %--------------------------%
    
    % -- Tx --
    [o_tx_s, cfg_s.tx_s] = m_qam_tx(cfg_s.tx_s);
    
    % -- Noise --
    n_sym = length(o_tx_s.sym_v);
    snr_db = cfg_s.ebno_db + 10*log10(log2(cfg_s.tx_s.M));
    pw_noise = var(o_tx_s.sym_v) / 10^(snr_db/10);
    ch_v = o_tx_s.sym_v + sqrt(pw_noise/2)*(randn(n_sym,1)+1j*randn(n_sym,1));
    
    % -- LW --
    fs = cfg_s.tx_s.BR;
    lw_sigma = sqrt(2*pi*cfg_s.lw*1/fs);
    lw_v = cumsum(lw_sigma*randn(n_sym,1));
    ch_v = ch_v .* exp(1j*lw_v);
    % -- BPS --
    i_bps_s.signal_v = ch_v.'; % El modulo recibe vectores fila
    o_bps_s = m_bps(i_bps_s, cfg_s.bps_s);
    
    rx_bps_v = o_bps_s.bps_output;
    bps_phase=o_bps_s.bps_phase;
    
    % CS corrector (BPS no corrige Cycle Slip)
    delay = (cfg_s.bps_s.filter_length-1)/2;
    rx_bps_v = rx_bps_v(1+delay:end);
    ak_v = o_tx_s.sym_v(1:length(rx_bps_v));
    
    if cfg_s.en_cs_corrector
    
        i_cs_cor_s.tx_signal_v = ak_v;
        i_cs_cor_s.rx_signal_v = rx_bps_v;
        o_cs_cor_s = m_cs_corrector(i_cs_cor_s, cfg_s.cs_cor_s);

        rx_bps_cs_v = o_cs_cor_s.signal_v;
        
    else
        rx_bps_cs_v = rx_bps_v;
    end
    
    % -- BER estimation for BPS --
    % Slicer
    ak_hat_bps_v = my_qam_slicer(rx_bps_cs_v, cfg_s.tx_s.M);
    
    % QAM to bits
    ak_bit_v = qamdemod(ak_v, cfg_s.tx_s.M, 'OutputType', 'bit');
    ak_hat_bit_v = qamdemod(ak_hat_bps_v,  cfg_s.tx_s.M, 'OutputType', 'bit');
    
    % Align 
    delay = finddelay(ak_bit_v, ak_hat_bit_v);
    if delay<0
        delay = abs(delay);
        ak_bit_v = ak_bit_v(1+delay:delay+length(ak_hat_bit_v));
    else
        ak_hat_bit_v = ak_hat_bit_v(1+delay:end);
        ak_bit_v = ak_bit_v(1:length(ak_hat_bit_v));
    end
    
    % BER estimation
    n_bit_errors = sum(ak_bit_v ~= ak_hat_bit_v);
    ber_est_bps = n_bit_errors / length(ak_bit_v);
    
    if n_bit_errors<100
        warning('BPS BER estimation may be not accurate')
    end
    
    % -- BER theo --
    ebno = 10^(snr_db/10) / log2(cfg_s.tx_s.M);
    ebno_db = 10*log10(ebno);
    ber_theo = berawgn(ebno_db, 'QAM', cfg_s.tx_s.M);
      % -- Calcule SNR de la se�al de entrada al bloque BPS --
signal_power_tx = mean(abs(o_tx_s.sym_v).^2);
noise_power_tx = var(ch_v ) - signal_power_tx;
snr_db_tx = 10 * log10(signal_power_tx / noise_power_tx);

% -- Calcule SNR de la se�al de salida del bloque BPS --
signal_power_rx = mean(abs(ak_hat_bps_v).^2);
noise_power_rx = var(rx_bps_cs_v) - signal_power_rx;
snr_db_rx = 10 * log10(signal_power_rx / noise_power_rx);

% -- Calcule la p�rdida de SNR --
snr_loss = snr(abs(o_tx_s.sym_v))- snr(abs(ak_hat_bps_v));
    %--------------------------%
    %          OUTPUT
    %--------------------------%
    
    if cfg_s.plots
        figure;
        plot(lw_v, '-', 'DisplayName','Ruido de Fase', 'LineWidth', 2);
        hold on;
        grid on;
        plot(unwrap(-bps_phase), '-', 'DisplayName','BSP', 'LineWidth',2);
        titulo = sprintf('Seguimiento del BSP \n Numeros de Fases: %d, Ruido de Fase (Khz): %.2f, Longitud del Filtro: %d',cfg_s.bps_s.n_phases, cfg_s.lw/1e3, cfg_s.bps_s.filter_length);
        hTitle = title(titulo);
        set(hTitle, 'FontSize', 14);
        ylabel('Amplitud');
        xlabel('Muestras');
        legend('Location', 'Best');
    end
    
    o_data_s.ber_theo = ber_theo;
    o_data_s.ber_est_bps = ber_est_bps;
    o_data_s.n_bit_errors = n_bit_errors;
    o_data_s.ebno_db = ebno_db;
    o_data_s.snr_loss = snr_loss;
    
end
    