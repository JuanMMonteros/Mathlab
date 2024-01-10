%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Santiago F. Leguizamon
% Created on   : July 2023
% Description  : Laser model
%-----------------------------------------------------------------------------%

function [o_data_s, o_cfg_s] = m_cs_corrector(i_data_s, i_cfg_s)

    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%
    
    cfg_s.window_length = 50    ;

    % --- Private --- % 

    
    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%

    if nargin > 0
        cfg_s = overwrite_parameters(i_cfg_s, cfg_s);
    end
    
    %--------------------------%
    %          ERRORS
    %--------------------------%
    
    %--------------------------%
    %     LOCAL VARIABLES
    %--------------------------%

    data_length = length(i_data_s.rx_signal_v);
    n_blocks = fix(data_length/cfg_s.window_length);
    rx_cs_fixed_v = zeros(n_blocks * cfg_s.window_length, 1);
    cs_phase = zeros(n_blocks, 1);
    cs_count = 0;
    last_phase = 0;
        
    %--------------------------%
    %          PROCESS
    %--------------------------%

    for block_idx = 1 : n_blocks
        
        slice = (block_idx - 1) * cfg_s.window_length + 1 : block_idx * cfg_s.window_length;
        rx_block_in = i_data_s.rx_signal_v(slice);
        tx_block_in = i_data_s.tx_signal_v(slice);
        
        min_mse = inf;
        phase_ok = 0;
        
        % Testing all the phases and saving the one which produces the smallest
        % MSE
        for phase_test_value = [0, pi/2, -pi/2, pi]
            
            block_test = rx_block_in .* exp(1j * phase_test_value);
            mse = mean(abs(block_test - tx_block_in).^2);
            if mse < min_mse
                min_mse = mse;
                phase_ok = phase_test_value;
            end
            
        end

        if phase_ok ~= last_phase
            cs_count = cs_count + 1;
        end
        
        last_phase = phase_ok;
        cs_phase(block_idx) = phase_ok; 
        rx_cs_fixed_v(slice) = rx_block_in.*exp(1j*phase_ok);
        
    end
        
    %--------------------------%
    %          OUTPUT
    %--------------------------%
    
    o_data_s.signal_v = rx_cs_fixed_v;
    o_data_s.cs_count = cs_count;
    o_data_s.cs_phase = cs_phase;
    o_cfg_s = cfg_s;
    
    %--------------------------%
    %           PLOTS
    %--------------------------%
    

end



