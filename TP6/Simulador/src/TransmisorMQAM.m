%-------------------------------------------------------------------------%
%                                 FULGOR
% Programmer(s): Daniel A. Juarez
% Created on   : October 2023
% Description  : Transmitter Filter Tx
%-------------------------------------------------------------------------%

function [o_data_s] = TransmisorMQAM(i_config_s)
    
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%
    BR = 32e9;
    fs = BR;

    config.BR = 32e9;    % Baud rate
    config.M = 4;        % Niveles de la modulacion
    config.L = 10e3;     % Largo de simu
    config.fs = BR;      % Frec de sampling
    config.T = 1/BR;     % Intervalo de tiempo entre dos simbolos consecutivos
    config.Ts = 1/fs;    % Tiempo entre 2 samples consecutivos a la salida del Tx
    
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
    M = config.M;
    L = config.L;

    %--------------------------%
    %         PROCESS
    %--------------------------%

    % Two symbols generation (+1,-1) for QPSK
    xsymb = qammod(randi([0 M-1], L,1), M);

    %--------------------------%
    %         OUTPUT
    %--------------------------%
    o_data_s.xsymb = xsymb;    
   
end

