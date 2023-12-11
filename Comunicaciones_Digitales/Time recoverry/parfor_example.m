clear; clc; close all;

% Default configuration
config_s.M = 4;
config_s.BR = 10e9;

% Parameters to sweep

rolloff_v = [0.1, 0.5, 0.7];
n_rolloff = length(rolloff_v);

ntaps_v = [10 20];
n_ntaps = length(ntaps_v);

step_v = [1e-6, 2e-6, 4e-6, 8e-6];
n_step = length(step_v);

% Sweep

for idx1 = 1:n_rolloff
    
    config_s.rolloff = rolloff_v(idx1);
    
    for idx2 = 1:n_ntaps
        
        config_s.ntaps = ntaps_v(idx2);
        
        parfor idx3 = 1:n_step
            
            cfg_s = config_s;
            cfg_s.step = step_v(idx3);
            
            function_print(cfg_s);
            
        end
        fprintf("\n");
    end
end

% Funcion para imprimir
function function_print(i_struct)

    fprintf(" - Rolloff = %.1f. Taps = %d. Step = %.1e\n", ...
                    i_struct.rolloff, i_struct.ntaps, i_struct.step);
    pause(3); % Pongo esta pausa para mostrar bien que se ejecuta a la vez

end
            
           