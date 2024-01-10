%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero 
% Created on   : September 2023
% Description  : CPR simulator
%-----------------------------------------------------------------------------%

function o_data_s = cpr_sim(i_cfg_s)

    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%
    
    % -- Tx --
    cfg_s.tx_s.n_sym                  = 1e5     ;   % Number of symbols
    cfg_s.tx_s.n_pil                  = 100     ;   % Number of pilots symbols
    cfg_s.tx_s.n_payload              = 127     ;   % Symbols between pilots
    
    cfg_s.tx_s.M                      = 16      ;   % Modulation levels
    cfg_s.tx_s.BR                     = 240e9   ;   % Baud rate [Bd]
   
    % -- Laser --
    cfg_s.lw                          = 0e3     ;   % Tx LW [Hz]
    
    % -- Noise --
    cfg_s.ebno_db                     = 10      ;   % EbNo [dB]
    
    % -- CPR --
    cfg_s.n_taps_cpr                  = 65      ;   % Taps to filter CPR error
    
    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%

    if nargin > 0
        cfg_s = overwrite_parameters(i_cfg_s, cfg_s);
    end

    clear i_cfg_s;
    
    %--------------------------%
    %          PROCESS
    %--------------------------%
    
    % -- Tx --
    [o_tx_s, cfg_s.tx_s] = m_qam_tx(cfg_s.tx_s);
    
    % -- Noise --
    n_sym = length(o_tx_s.sym_v);
    snr_db = cfg_s.ebno_db + 10*log10(log2(cfg_s.tx_s.M));
    pw_noise = var(o_tx_s.sym_v) / 10^(snr_db/10);
    rx_v = o_tx_s.sym_v + sqrt(pw_noise/2)*(randn(n_sym,1)+1j*randn(n_sym,1));
    
    % -- LW --
    fs = cfg_s.tx_s.BR;
    lw_sigma = sqrt(2*pi*cfg_s.lw*1/fs);
    lw_v = cumsum(lw_sigma*randn(n_sym,1));
    rx_v = rx_v .* exp(1j*lw_v);
    
    % -- CPR --
    ak_v = o_tx_s.sym_v;
    n_ak = length(ak_v);
    n_rx = length(rx_v);
    rx_v = rx_v(n_rx-n_ak+1:end);
    fifo_len = ceil((cfg_s.tx_s.n_payload+1)/2) + ...
                        floor(cfg_s.n_taps_cpr/2) * (cfg_s.tx_s.n_payload+1);

    if cfg_s.tx_s.n_pil>0
    
        % Align to start with a pilot
        pil_ext_v = upsample(o_tx_s.pilot_v, cfg_s.tx_s.n_payload+1);
        ak_cor_v = ak_v(1:2*length(pil_ext_v)-1);
        delay = finddelay(pil_ext_v, ak_cor_v);
        ak_v = ak_v(delay+1:end);
        rx_v = rx_v(delay+1:end);
        n_ak = length(ak_v);

        % Variables
        h_cpr_v = 1/cfg_s.n_taps_cpr * ones(1, cfg_s.n_taps_cpr);
        h_zoh_v = ones(cfg_s.tx_s.n_payload+1, 1);
        cpr_error_v = zeros(ceil(n_ak/cfg_s.tx_s.n_payload),1);

        % Loop
        for idx = 1:n_ak

            if mod(idx-1, cfg_s.tx_s.n_payload+1)==0

                angle_dif = angle(rx_v(idx)) - angle(ak_v(idx));
                if angle_dif > pi
                    angle_dif = angle_dif - 2*pi;
                elseif angle_dif < - pi
                    angle_dif = angle_dif + 2*pi;
                end

                cpr_error_v(floor(idx/(cfg_s.tx_s.n_payload+1)+1)) = angle_dif;

            end
        end
        cpr_error_fil_v = filter(h_cpr_v, 1, unwrap(cpr_error_v));
        cpr_phase_v = filter(h_zoh_v, 1, ...
                            upsample(cpr_error_fil_v, cfg_s.tx_s.n_payload+1));

        cpr_phase_align_v = cpr_phase_v(fifo_len+1:end);
        if length(cpr_phase_align_v) < n_ak
            n_ak = length(cpr_phase_align_v);
        else
            cpr_phase_align_v = cpr_phase_align_v(1:n_ak);
        end
        rx_v = rx_v(1:n_ak) .* exp(-1j*cpr_phase_align_v);
        ak_v = ak_v(1:n_ak);

        % Take out pilots to count BER
        rx_v(1:cfg_s.tx_s.n_payload+1:end) = [];
        ak_v(1:cfg_s.tx_s.n_payload+1:end) = [];
        
    end
    
    % -- BER estimation --
    % Slicer
    ak_hat_v = my_qam_slicer(rx_v, cfg_s.tx_s.M);
    
    % QAM to bits
    ak_bit_v = qamdemod(ak_v, cfg_s.tx_s.M, 'OutputType', 'bit');
    ak_hat_bit_v = qamdemod(ak_hat_v,  cfg_s.tx_s.M, 'OutputType', 'bit');
    
    % BER estimation
    n_bit_errors = sum(ak_bit_v ~= ak_hat_bit_v);
    ber_est = n_bit_errors / length(ak_bit_v);
    
    if n_bit_errors<100
        warning('BER estimation may be not accurate')
    end
    
    % -- BER theo --
    ebno = 10^(snr_db/10) / log2(cfg_s.tx_s.M);
    ebno_db = 10*log10(ebno);
    ber_theo = berawgn(ebno_db, 'QAM', cfg_s.tx_s.M);
    
    %--------------------------%
    %          OUTPUT
    %--------------------------%
    
    o_data_s.ber_theo = ber_theo;
    o_data_s.ber_est = ber_est;
    o_data_s.n_bit_errors = n_bit_errors;
    o_data_s.ebno_db = ebno_db;
    
end
    