%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero
% Created on   : July 2023
% Description  : QAM Tx
%-----------------------------------------------------------------------------%

function [o_data_s, o_cfg_s] = m_qam_tx(i_cfg_s)

    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    cfg_s.n_sym            = 1e6;       % Number of symbols
    cfg_s.n_pil            = 100;       % Number of pilots symbols
    cfg_s.n_payload        = 127;       % Symbols between pilots
    
    cfg_s.M                = 64;        % Modulation levels
    cfg_s.BR               = 240e9;     % Baud rate
    
    cfg_s.OVR              = 2;         % Oversampling rate
    cfg_s.rolloff          = 0.5;       % PS rolloff
    cfg_s.n_taps           = 31;        % PS number of taps
    cfg_s.h_ps_v           = NaN;       % PS taps
    cfg_s.ps_selec         = 1;         % PS sel. 0:OFF;1:RRC;2:RC;3:PE

    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%

    if nargin > 0
        cfg_s = overwrite_parameters(i_cfg_s, cfg_s);
    end

    clear i_cfg_s;
    
    %--------------------------%
    %           ERRORS
    %--------------------------%
    
    if mod(log2(cfg_s.M), 2) ~= 0 && cfg_s.M~=2
        error('The modulation levels must be M=2^k with even k (or k=1)')
    end
    
    if mod(cfg_s.n_taps, 2) == 0
        error('PS number of taps must be odd')
    end
    
    %--------------------------%
    %   CONSTANTS & VARIABLES
    %--------------------------%
    
    % Oversampling
    [OVR_num, OVR_den] = rat(cfg_s.OVR);
    
    % Pulse shaping
    ps_span = (cfg_s.n_taps-1)/2;
    if cfg_s.ps_selec == 0
        h_ps_v = zeros(cfg_s.n_taps,1);
        h_ps_v(ps_span) = 1;
    elseif cfg_s.ps_selec == 1
        h_ps_v = root_raised_cosine(cfg_s.BR/2, OVR_num*cfg_s.BR, ...
                                            cfg_s.rolloff, cfg_s.n_taps, 0);
    elseif cfg_s.ps_selec == 2
        h_ps_v = raised_cosine(cfg_s.BR/2, OVR_num*cfg_s.BR, ...
                                            cfg_s.rolloff, cfg_s.n_taps, 0);
    elseif cfg_s.ps_selec == 3
        h_ps_v = cfg_s.h_ps_v;
    end
    
    h_ps_v = OVR_num * h_ps_v / sum(h_ps_v);
    h_grpdelay_v = grpdelay(h_ps_v);
    ps_delay = round(h_grpdelay_v(1));
    
    %--------------------------%
    %          PROCESS
    %--------------------------%
    
    % Symbols generation
    randi_v = randi([0 cfg_s.M - 1], cfg_s.n_sym, 1);
    sym_payload_v = qammod(randi_v, cfg_s.M);

    % Pilots generation (QPSK)
    if cfg_s.n_pil == 0
        
        pilot_v = [];
        sym_v = sym_payload_v;
        
    else
        
        pilot_v = (sqrt(cfg_s.M)-1) * qammod(randi([0 3], cfg_s.n_pil, 1), 4);
        
        % Pilots insertion
        n_total_pilots = ceil(cfg_s.n_sym/cfg_s.n_payload);
        sym_v = zeros(cfg_s.n_sym + n_total_pilots, 1); 
        idx_pil = 1; n_block = cfg_s.n_payload;

        for idx = 0 : n_total_pilots-2

            sym_v(idx*(n_block+1) + 1) = pilot_v(idx_pil);
            sym_v(idx*(n_block+1) + 2 : idx*(n_block+1) + n_block + 1) = ...
                            sym_payload_v(idx*n_block + 1 : idx*n_block + n_block);

            idx_pil = idx_pil + 1;
            if idx_pil > cfg_s.n_pil
                idx_pil = 1;
            end

        end

        if sym_v(end) == 0
            idx = idx + 1;
            sym_v(idx*(n_block+1) + 1) = pilot_v(idx_pil);
            sym_v(idx*(n_block+1) + 2 : end) = ...
                                            sym_payload_v(idx*n_block + 1 : end);
        end
        
    end
        
    % Power normalization
    rx_mf_gain = max(real(sym_v));
    sym_norm_v = sym_v / rx_mf_gain;
    
    % Upsampling 
    sym_up_v = upsample(sym_norm_v, OVR_num);
    
    % Filter
    sym_fil_v = filter(h_ps_v, 1, [sym_up_v; zeros(ps_delay, 1)]);
    sym_fil_v = sym_fil_v(1+ps_delay:end);
    
    % Downsampling 
    tx_v = downsample(sym_fil_v, OVR_den);
    
    % CMA rate
    if cfg_s.n_sym < 1000
        cma_ref = sqrt(mean(abs(sym_v(1:cfg_s.n_sym)).^4) / ...
                                    mean(abs(sym_v(1:cfg_s.n_sym)).^2));
    else
        cma_ref = sqrt(mean(abs(sym_v(1:1000)).^4) / ...
                                    mean(abs(sym_v(1:1000)).^2));
    end
    
    %--------------------------%
    %          OUTPUT
    %--------------------------%

     o_data_s.tx_v = tx_v;
    
    o_data_s.sym_v = sym_v;
    o_data_s.pilot_v = pilot_v;
    
    o_data_s.h_ps_v = h_ps_v;
    
    o_data_s.cma_ref = cma_ref;
    o_data_s.rx_mf_gain = rx_mf_gain;
    
    o_cfg_s = cfg_s;

    %--------------------------%
    %           PLOTS
    %--------------------------%

end
