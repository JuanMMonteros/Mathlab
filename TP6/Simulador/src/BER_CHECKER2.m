%-------------------------------------------------------------------------%
%                                   FULGOR
% Programmer(s): Daniel A. Juarez
% Created on   : October 2023
% Description  : BER checker
%-------------------------------------------------------------------------%

function [o_data_s] = BER_CHECKER2(i_ber_s)
    
    %--------------------------%
    %   CONSTANTS & VARIABLES
    %--------------------------%    
    M = i_ber_s.M;                   % Orden de modulación
    EbNo_dB = i_ber_s.EbNo_dB;
    xsymb = i_ber_s.xsymb;           % Constelación salida Tx
    ak_hat = i_ber_s.ak_hat;         % Simbolos recuperados   
    pll_output = i_ber_s.pll_output; % Salida del PLL        
    Ldata = i_ber_s.Ldata;           % Longitud de los datos

    %--------------------------%
    %          PROCESS
    %--------------------------%
    % 1. Alinear
    d0 = finddelay(xsymb, pll_output);
    guard0 = fix(Ldata*.25);
    guard1 = 1e3;
    a_hat_align = ak_hat (1+d0+guard0:end-guard1);
    pll_output_align = pll_output (1+d0+guard0:end-guard1);
    x_align = xsymb (1+guard0:end-guard1-d0);
    
    % 2. OPCION 1 (millenial, facil)
    simbolos_totales = length(x_align);
    simbolos_errados = sum(x_align ~= a_hat_align);
    symbol_error_rate = simbolos_errados / simbolos_totales;
    aprox_ber_sin_corregir_cs = 1/log2(M) * symbol_error_rate;
    
    % %%
    % bits_tx = demodulador_M_QAM(x_align, M);
    % bits_rx = demodulador_M_QAM(a_hat_align, M);
    % [errors, ber] = biterr(bits_tx, bits_rx);
    % ber;
    
    
    % 
    % %% Correccion dinamica de los cycle slips
    % Hacerlo con los simbolos alineados
    % orx son los simbolos a la salida del BPS (sin slicer)
    % otx son los simbolos transmitidos
    % orx y otx estan alineados
    orx = pll_output_align;
    otx = x_align;
    WINDOW_LEN = 500;
    Ldata = length(orx);
    nblocks = fix(Ldata/WINDOW_LEN);
    
    orx_cs_fixed = zeros(nblocks*WINDOW_LEN,1);
    cs_phase = zeros(nblocks,1);
    cs_count =0;
    last_phase=0;
    mse_log = zeros(nblocks, 4);
    phase_test_list = [0, pi/2, -pi/2, pi];
    
    for nblock=1:nblocks
        slice = (nblock-1)*WINDOW_LEN+1: nblock*WINDOW_LEN;
        rx_block_in = orx(slice);
        tx_block_in = otx(slice);
        
        min_mse = inf;
        phase_ok = 0;
        for np = 1:4
            phase_test = phase_test_list(np);
            block_test = rx_block_in.*exp(1j*phase_test);
            mse = mean(abs(block_test-tx_block_in).^2);
            mse_log(nblock, np) = mse;
            if mse < min_mse
                min_mse = mse;
                phase_ok=phase_test;
            end
        end
        if phase_ok ~= last_phase
            cs_count = cs_count +1;
        end
        last_phase = phase_ok;
        cs_phase(nblock) = phase_ok; 
        orx_cs_fixed(slice) = rx_block_in.*exp(1j*phase_ok);
    end
    
    figure
    plot(cs_phase/(pi/2))
    ylim([-2,2])
    
    Ldata = length(orx_cs_fixed);
    orx_cs_fixed_slicer = zeros(Ldata,1);
    for m=1:Ldata
        orx_cs_fixed_slicer(m) = slicer(orx_cs_fixed(m), M);
    end
    
    simbolos_totales = length(orx_cs_fixed_slicer);
    simbolos_errados = sum(x_align ~= orx_cs_fixed_slicer);
    symbol_error_rate = simbolos_errados / simbolos_totales;
    aprox_ber_corrigiendo_cs = 1/log2(M) * symbol_error_rate;
    
    ber_teo = berawgn(EbNo_dB, 'qam', M);

    %--------------------------%
    %         OUTPUT
    %--------------------------%  
    o_data_s.aprox_ber_sin_corregir_cs = aprox_ber_sin_corregir_cs;
    o_data_s.aprox_ber_corrigiendo_cs = aprox_ber_corrigiendo_cs;
    o_data_s.ber_teo = ber_teo;
end
