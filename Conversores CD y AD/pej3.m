clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

fs = 400e6;
Ts = 1/fs;
OVR = 256;
fs_ch = fs * OVR;
fz = 15;    % Plots Font Size

n_samples = 1e2;
% Sim length in digital domain
omega_0 = 2*pi*50e6;   % Tone frequency
T_tone = 2*pi/(omega_0); % Periodo del Tono 

% Time vector
t_up_v = (0:n_samples*OVR-1) * 1/fs_ch;
t_v = (0:n_samples-1) * Ts;

% Ideal signal in digital domain
tone_dig = 2 * cos(omega_0 * t_v);

% Signal in continuous domain
tone_cont = 2 * cos(omega_0 * t_up_v);  

%% DELAY ANALÓGICO
% Se pretende atrasar 40 muestras la señal
delay_analogico = [zeros(1, 40), tone_cont(1:end-40)];
% -- Downsampling --
tone_dig_CD = delay_analogico(1:OVR:end);

%% DELAY DIGITAL 
% Se pretende retrasar la señal una fracción de muestra
retardo = 390.625e-12;  % En segundos
delay = retardo/Ts;     % Esto representa una fracción de muestra
% Genero el filtro de la eq 4.65 del libro
ntaps = 31;             % Mantenerlo impar para simplificar
group_delay = (ntaps-1)/2; % Este delay hay que compensarlo porque representa la cola de la convolución
nline = -(ntaps-1)/2:(ntaps-1)/2;
h = sinc(nline - delay);
tone_dig1 = transpose(tone_dig);
tone_dig1 = [tone_dig1; zeros(group_delay, 1)];
tone_dig_delay = filter(h, 1, tone_dig1); % Agrego zeros para mantener el largo de la señal filtrada
tone_dig_delay = tone_dig_delay(group_delay+1:end); % Corrijo el retardo de grupo

%% PLOT EN TIEMPO
figure
% Subplot de arriba: se muestra la señal discreta y la señal continua antes del delay
subplot(2, 1, 1);
p1 = plot(t_up_v/1e-6, tone_cont ,'-r','Linewidth',1);
p1.MarkerFaceColor = p1.Color;
p1.MarkerEdgeColor = p1.Color;
hold on;
p2 = plot(t_v/1e-6, tone_dig ,'ob','Linewidth',1.5);
xlabel('Tiempo (us)');
ylabel('Amplitud');
title('Señales discretas y continuas antes del delay');
legend([p1, p2], 'Señal continua', 'Señal discreta');
xlim([T_tone/1e-6, 3*T_tone/1e-6]); 

% Subplot de abajo: se muestra la señal discreta y la señal continua después del delay
subplot(2, 1, 2);
xlabel('Tiempo (us)');
ylabel('Amplitud');
p3 = plot(t_up_v/1e-6, delay_analogico ,'--k','Linewidth',1);
p3.MarkerFaceColor = p3.Color;
p3.MarkerEdgeColor = p3.Color;
hold on;
p4 = plot(t_v/1e-6, tone_dig_delay ,'ob','Linewidth',1.5);
title('Señales discretas después del delay (usando sinc)');
legend([p3, p4], 'Señal continua después del delay', 'Señal discreta después del delay');
% Mostrar solo dos periodos para mayor claridad
xlim([T_tone/1e-6, 3*T_tone/1e-6]);
% Ajustar los subplots para una mejor visualización
sgtitle('Simulación de Tono discreto y Tono continuo con delay fraccional usando filtro sinc');

%% FFT variables
NFFT = 16 * n_samples;
f_v = (-NFFT/2:NFFT/2-1) * fs/NFFT;
f_up_v = (-NFFT*OVR/2:NFFT*OVR/2-1) * fs_ch/(NFFT*OVR);

% FFTs
X_d = fftshift(abs(fft(tone_dig, NFFT))/n_samples);
X_c = fftshift(abs(fft(tone_cont, NFFT*OVR))/(n_samples*OVR));
X_d_delay = fftshift(abs(fft(tone_dig_delay, NFFT))/(n_samples));
X_c_delay = fftshift(abs(fft(delay_analogico, NFFT*OVR))/(n_samples*OVR));
H_v = fftshift(abs(fft(h, NFFT*OVR)));

% Respuesta en frecuencia de la señal analógica
X_analog = fftshift(abs(fft(delay_analogico, NFFT*OVR))/(n_samples*OVR));

% Respuesta en frecuencia de la señal digital antes y después del remuestreo
X_d_digital = fftshift(abs(fft(tone_dig, NFFT*OVR))/(n_samples*OVR));
X_d_digital_delay = fftshift(abs(fft(tone_dig_delay, NFFT))/(n_samples));

% Freq
figure
plot(f_up_v/1e6, X_analog,'-r','Linewidth',2);
hold on;
plot(f_v, X_d ,'--b','Linewidth',2);
plot(f_v, X_d_delay ,'-.c','Linewidth',1);
plot(f_up_v, H_v ,'-.m','Linewidth',1);
tit = sprintf('Respuesta en Frecuencia');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Frecuencia (MHz)', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitud', 'Interpreter','latex','FontSize', fz);
legend({'Señal analógica','Señal digital antes del remuestreo','Señal digital después del remuestreo','Respuesta del filtro'}, 'Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [550 50 800 500],'Color', 'w');
