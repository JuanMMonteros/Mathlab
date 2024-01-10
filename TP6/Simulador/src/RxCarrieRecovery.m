%-------------------------------------------------------------------------%
%                                 FULGOR
% Programmer(s): Daniel A. Juarez
% Created on   : October 2023
% Description  : Reception FSE (Fractionaly Spacer equalizer)
%-------------------------------------------------------------------------%

function [o_data_s] = RxCarrieRecovery(i_rx_s, i_config_s)
    
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    config.Kp = 15e-2;      % Constante Kp
    config.Ki = 15e-2/1000; % Constante Ki 
    config.M = 4;           % Nivel de modulación

    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%

    fn = fieldnames(i_config_s);
    for k = 1:numel(fn)
        if isfield(config,(fn{k}))==1
            config.(fn{k})= i_config_s.(fn{k});
        else
            error("%s: Parametro del simulador no valido", fn{k})
        end
    end

    %--------------------------%
    %         VARIABLES
    %--------------------------%
    % -- Variables de entrada -- %
    symb_in = i_rx_s.ch_out_symb;
    
    Kp = config.Kp;       
    Ki = config.Ki; 
    M = config.M;           

    %--------------------------%
    %         PROCESS
    %--------------------------%

    % -- PLL dirigido por decisiones -- %
    Ldata = length(symb_in);
    phase_error = zeros(Ldata,1);
    nco_output = zeros(Ldata,1);
    integral_branch = zeros(Ldata,1);
    pll_output = zeros(Ldata,1);
    
    for m=2:Ldata
        xpll = symb_in(m);
        derot_x = xpll*exp(-1j*nco_output(m-1));
        pll_output(m) = derot_x;
        
        a_hat = slicer(derot_x, M);
        phase_error(m) = angle(derot_x .* conj(a_hat));
        prop_error = phase_error(m) * Kp;
        integral_branch(m) = integral_branch(m-1) + Ki*phase_error(m);
        nco_output(m) = nco_output(m-1) + prop_error +integral_branch(m);
    end

    % -- Slicer -- %
    z = pll_output;
    Ldata = length(z);
    ak_hat = zeros(Ldata,1);
    for m=1:Ldata
        ak_hat(m) = slicer(z(m), M);
    end
    
    %--------------------------%
    %         OUTPUT
    %--------------------------%
    o_data_s.pll_output = pll_output;           % Salida del PLL
    o_data_s.nco_output = nco_output;           % Salida del NCO
    o_data_s.integral_branch = integral_branch; % Salida de la rama integral
    o_data_s.ak_hat = ak_hat;                   % Coef ak 
    o_data_s.phase_error = phase_error;         % Error de fase
    o_data_s.Ldata = Ldata;                     % Longitud del dato
    
end
