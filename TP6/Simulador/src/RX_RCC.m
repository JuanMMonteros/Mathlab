function [o_data_s] = RX_RCC(i_config_s,rx_data)

 %--------------------------%
 %     DEFAULT SETTINGS
 %--------------------------%

    config.Kp = 0.05;      
    config.Ki = config.Kp/1000;  
    config.M = 4;           

    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%

    fn = fieldnames(i_config_s);
    for k = 1:numel(fn)
        if isfield(config,(fn{k}))==1
            config.(fn{k})= i_config_s.(fn{k});
        %else
           % error("%s: Parametro del simulador no valido", fn{k})
        end
    end

    %--------------------------%
    %         VARIABLES
    %--------------------------%
    % -- Variables de entrada -- %   
    Kp = config.Kp;       
    Ki = config.Ki; 
    M = config.M;           

    %--------------------------%
    %         PROCESS
    %--------------------------%
    
     Ldata = length(rx_data);
    phase_error = zeros(Ldata,1);
    nco_output = zeros(Ldata,1);
    integral_branch = zeros(Ldata,1);
    pll_output = zeros(Ldata,1);
    
    for m=2:Ldata
        xpll = rx_data(m);
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