function [o_data_s] = AGC(i_config_s,rx)
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    config.target = 0.3; 
    config.adc_phase= 1;
    config.NOS= 4;
    config.N=2;

    
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
    adc_phase=config.adc_phase;
    target= config.target;
    N=config.N;
    NOS=config.NOS;
    
    %--------------------------%
    %         PROCESS
    %--------------------------%
    rx = rx(adc_phase: NOS/N : end);
    metric = std(rx);
    agc_gain = target/metric;
    rx_norm = rx.*agc_gain;
    o_data_s.rx_norm=rx_norm;
    
   