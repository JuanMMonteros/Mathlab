%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        TP6 EJERCICIO 5                          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
clear all 
clc 

%%  -- Modelado de un PLL Tipo 2 --

% -- Parametros de simulación -- %
BR = 16e9;      % Baurate 
T = 1/BR;       % Periodo de muestreo
fs = BR;        % Frec de samplin
L = 50e3;       % Longitud de simulacion
Lat = 50;     % Latencia de simulación

% -- Valores de Kp y Ki a barrer -- %
Kp = 0.02;   
Ki = Kp/1000;

% -- Parametros para los plots -- %
t = [0:L-1].*T; % base de tiempo para los plots  
fz = 15;
position_and_size = [100 100 800 600];
saves_plots = 1; % Enable save plots 1 -> ON, 0 -> OFF

%%
% -- Path para salvar los plots -- %

figur_dir = mfilename('fullpath');    % Genera directorio de salida
figur_dir = figur_dir(1:end-length(mfilename));
figur_dir = [figur_dir, 'Plots/'];

if ~exist(figur_dir,'dir')
    mkdir(figur_dir);
end


%% Punto a) 

% -- Step -- %
delay_init = 1000;
tita0 = 20; % Grados
tita_in = [zeros(delay_init,1); tita0/180*pi*ones(L-delay_init,1)]; % Escalon de fase

% -- Variables de logueo -- %
error_input = zeros(L, 1);
error_prop = zeros(L, 1);
error_integral = zeros(L, 1); 
NCO_in = zeros(L, 1);
NCO_out = zeros(L, 1);
tital_delay = zeros(L, 1);
DELAY_LINE = zeros(Lat+1,1);

for n=2:L
    % Error de la entrada
    error_input(n) = tita_in(n) - tital_delay(n-1);
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
    tital_delay(n) = DELAY_LINE(1);
end

% -- Plots -- %
% Escalon de entrada vs salida
figure
plot(t, tita_in, '--b', 'Linewidth',2);
hold all
plot(t, NCO_out, '-r', 'Linewidth',2);
grid on
xlim([0, 0.2e-6]);
tit = sprintf('Respuesta al Escalon');
saveName = sprintf('Ejerc5_a1)_StepFrecResp.jpg');
name1 = sprintf('Entrada PLL');
name2 = sprintf('Salida PLL');
legend(name1,name2,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
title(tit,'Interpreter','latex', 'FontSize', fz);
xlabel('Time [s]')
ylabel('Amplitudes [rad]')
set(gcf, 'Position',position_and_size,'Color', 'w');
if saves_plots              
    saveas(gcf,[figur_dir, saveName]);          
end

% Error estacionario del sistema 
figure
plot(t, error_input, '-b', 'Linewidth',2)
hold all
grid on
xlim([0, 0.2e-6]);
tit = sprintf('System Error Evolution');
saveName = sprintf('Ejerc5_a2)_EssEvol.jpg');
name = sprintf('Ess evolución');
legend(name,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
title(tit,'Interpreter','latex', 'FontSize', fz);
xlabel('Time [s]')
ylabel('Amplitudes [rad]')
set(gcf, 'Position',position_and_size,'Color', 'w');
if saves_plots              
    saveas(gcf,[figur_dir, saveName]);          
end

%% Punto b) 

% -- Ramp -- %
f0 = 200e6;      % Frec en Hz
w0 = 2*pi*f0;    % Frec en Omega
t = [0:L-1].*T;  % Tiempo de simulación
tita_in = w0.*t; % Entrada Rampa        

% -- Variables de logueo -- %
error_input = zeros(L, 1);
error_prop = zeros(L, 1);
error_integral = zeros(L, 1); 
NCO_in = zeros(L, 1);
NCO_out = zeros(L, 1);
tital_delay = zeros(L, 1);
DELAY_LINE = zeros(Lat+1,1);

for n=2:L
    % Error de la entrada
    error_input(n) = tita_in(n) - tital_delay(n-1);
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
    tital_delay(n) = DELAY_LINE(1);
end

% -- Plots -- %
% Rampa de entrada vs salida
figure
plot(t, tita_in, '--b', 'Linewidth',2);
hold all
plot(t, NCO_out, '-r', 'Linewidth',2);
grid on
xlim([0, 0.5e-6]);
tit = sprintf('Respuesta a un escalón de entrada');
saveName = sprintf('Ejerc5_b1)_StepFrecResp.jpg');
name1 = sprintf('Step in');
name2 = sprintf('Step out');
legend(name1,name2,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
title(tit,'Interpreter','latex', 'FontSize', fz);
xlabel('Time [s]')
ylabel('Amplitudes [rad]')
set(gcf, 'Position',position_and_size,'Color', 'w');
if saves_plots              
    saveas(gcf,[figur_dir, saveName]);          
end

% Error estacionario del sistema 
figure
plot(t, error_input, '-b', 'Linewidth',2)
hold all
grid on
xlim([0, 0.5e-6]);
tit = sprintf('Evolución del error en el sistema');
saveName = sprintf('Ejerc5_b2)_EssEvol.jpg');
name = sprintf('Ess evolución');
legend(name,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
title(tit,'Interpreter','latex', 'FontSize', fz);
xlabel('Time [s]')
ylabel('Amplitudes [rad]')
set(gcf, 'Position',position_and_size,'Color', 'w');
if saves_plots              
    saveas(gcf,[figur_dir, saveName]);          
end 

% Rama integral vs Rama Proporcional 
figure
plot(t, error_prop, '-b', 'Linewidth',2)
hold all
plot(t, error_integral, '-r', 'Linewidth',2)
grid on
tit = sprintf('Error proporcional vs Integral');
saveName = sprintf('Ejerc5_b3)_ErrPropVsErrIntrg.jpg');
name1 = sprintf('Err prop.');
name2 = sprintf('Err intrg.');
legend(name1,name2,'Location', 'eastoutside','Interpreter','latex','FontSize', fz-5);
title(tit,'Interpreter','latex', 'FontSize', fz);
xlabel('Time [s]')
ylabel('Amplitudes [rad]')
set(gcf, 'Position',position_and_size,'Color', 'w');
if saves_plots              
    saveas(gcf,[figur_dir, saveName]);          
end 

%% Punto c)

% -- Variables de logueo -- %
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
semilogx(fn_values, pll_out,'-r', 'LineWidth', 2);
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



