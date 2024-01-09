
close all 
clear all 
clc 



BR = 32e9;     
T = 1/BR;     
fs = BR;        
L = 50e3;      
Lat = 50;    
Kp = 0.05;   
Ki = Kp/1000;
t = [0:L-1].*T; 


%% Generacion de la entrada
input = 'Step'
switch(input)
    case 'Step'
        delay_init = 1000;
        tita0 = 20; % Grados
        tita_in = [zeros(delay_init,1); tita0/180*pi*ones(L-delay_init,1)]; % Escalon de fase
    case 'Ramp'
        f0 = 200e6; %Hz
        w0 = 2*pi*f0;
        t = [0:L-1].*T;
        tita_in = w0.*t;
end

%% Modelado del PLL
error = zeros(L,1);
error_prop = zeros(L,1);
error_int = zeros(L,1);
nco_in = zeros(L,1);
tita_out = zeros(L,1);
tita_delay = zeros(L,1);

DELAY_LINE = zeros(Lat+1,1); % FIFO para modelar la latencia

for n=2:L
    error(n) = tita_in(n) - tita_delay(n-1);
    error_prop(n) = Kp*error(n);
    
    % Calculo del error integral
    % Poco eficiente
    % if n<=1
    %     error_int(n) = Ki*error(n);
    % else
    %     error_int(n) = error_int(n-1)+Ki*error(n);
    % end
    % Mas eficiente es arrancar el for desde 2
    error_int(n) = error_int(n-1)+Ki*error(n);
    nco_in(n) = error_prop(n)+error_int(n);
    tita_out(n) = tita_out(n-1)+nco_in(n);
    
    % Modelado del retardo
    % Asumo que los datos entran por el final y salen por el frente
    % Primero roto, luego cargo el nuevo valor al final
    DELAY_LINE = [DELAY_LINE(2:end); tita_out(n)]; 
    % Finalmente, logueo la senial retardada
    tita_delay(n) = DELAY_LINE(1); 
end


%% 
%Plots

t = [0:L-1].*T;

% Plot respuesta a la entrada
figure
plot(t, tita_in, '--b')
hold all
plot(t, tita_out, '-r', 'Linewidth',2);
xlabel('Time [s]')
ylabel('Amplitudes [rad]')
grid on
legend('Input', 'Output')
title('Response to input')

% Error estacionario
figure
plot(t, error, '-b', 'Linewidth',2)
hold all
xlabel('Time [s]')
ylabel('Amplitude [rad]')
grid on
title('System Error Evolution')


% Errores 
figure
plot(t, error_prop, '-b', 'Linewidth',2)
hold all
plot(t, error_int, '-r', 'Linewidth',2)
plot(t, nco_in, '-g', 'Linewidth',2)

xlabel('Time [s]')
ylabel('Amplitude [rad]')
grid on
legend('Error Prop.', 'Error Int.', 'NCO In')
title('Proporcional and Integral Error')

% Error Integral escalado 
figure
plot(t, error_int*BR/2/pi, '-r', 'Linewidth',2)
xlabel('Time [s]')
ylabel('Frequency [Hz]')
grid on
title('Carrier Offset in Integrator Branch')
%%
error_input = zeros(L, 1);
error_prop = zeros(L, 1);
error_integral = zeros(L, 1); 
NCO_in = zeros(L, 1);
NCO_out = zeros(L, 1);
tita_delay = zeros(L, 1);
DELAY_LINE = zeros(Lat+1,1);
frec =logspace(log10(1e6), log10(fs/2), 1024); % Barrido de frec.
pll_out = zeros(length(frec),1); % Salida del sistema 
t = [0:L-1].*T; % Tiempo de simulacion

n_freq_pos = 2^16;  
n_freq_frac = 2^16;
n_v = [0, logspace(-log10(n_freq_frac),log10(n_freq_pos),n_freq_pos+n_freq_frac-1)];
f_v = n_v * BR/2/n_freq_pos;
wd_v = n_v * pi/n_freq_pos;

G_th_v = (Kp+Ki.*1./(1-exp(-1j*wd_v))) * 1./(1-exp(-1j.*wd_v));
H_th_v = G_th_v ./ (1+G_th_v.*exp(-1j*Lat*wd_v));
H_th_db_v = 20*log10(abs(H_th_v));

% -- Bode estimado -- %
for i = 1:length(frec)
    % -- Barrido de senoidales -- %
    fn = frec(i);
    wn = 2*pi*fn;
    tita0 = 15;
    t = [0:L-1].*T;
    tita_in = tita0/180*pi * sin(wn.*t);
    
    % -- Modelado del pll -- %
    for n=2:L
        % Error de la entrada
        error_input(n) = tita_in(n) - tita_delay(n-1);
        % Error de Kp
        error_prop(n) = Kp*error_input(n);
        % Error en Ki
        error_integral(n) = error_integral(n-1)+Ki*error_input(n);
        % Entrada del NCO
        NCO_in(n) = error_prop(n) + error_integral(n);
        % Salida del NCO
        NCO_out(n) = NCO_out(n-1) + NCO_in(n);
        % Modelado del retardo 
        DELAY_LINE = [DELAY_LINE(2:end); NCO_out(n)];
        tita_delay(n) = DELAY_LINE(1);
    end
    % Se obtiene los valores de salida
    pll_out(i) = 10*log10(var(NCO_out)/var(tita_in));

end
%%
% -- Plot Bode Teorico vs estimado -- % 
figure
semilogx(frec, pll_out,'-r', 'LineWidth', 2);
hold all 
semilogx(f_v, H_th_db_v, '--b', 'LineWidth', 2);
grid on 
xlabel('Frecuencia [Hz]');
ylabel('Amplitude [dB]');
xlim([1e6 1e9]);
title('Bode Teorico vs Estimado');
xlabel('Time [s]')
ylabel('Amplitudes [rad]')
set(gcf, 'Position',position_and_size,'Color', 'w');