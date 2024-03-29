
clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window
%% Parameters
fs=400e6;
Ts=1/fs;
OVR=256;
fs_ch=fs*OVR;
fz = 15;                % Plots Font Size
n_samples = 1e2;
% Sim length in digital domain
omega_0 = 2*pi*150e6;   % Tone frequency
T_tone=2*pi/(omega_0);% Periodo del Tono 

% Time vector
t_up_v = (0:n_samples*OVR-1)*1/fs_ch;
t_v = (0:n_samples-1)*Ts;

% Ideal signal in digital domain
tone_dig = 2*cos(omega_0*t_v); 

% Signal in contiuous domain
tone_cont = 2*cos(omega_0*t_up_v);  

%% DElAY ANALOGICO
% se prentede atrazar 40 muestras la senlas  
delay_analogico = [zeros(1, 40), tone_cont(1:end-40)];
% -- Downsampling --
tone_dig_CD = delay_analogico(1:OVR:end);   



%% DELAY DIGITAL 
% Se pretende retrasar la senial una fraccion de muestra
retardo=390.625e-12;
delay=retardo/Ts; % esto representa una fraccion de muestra
% Genero el filtro de la eq 4.65 del libro
ntaps=31; % Mantenerlo impar para simplificar
group_delay = (ntaps-1)/2; % Este delay hay que compensarlo porque representa la cola de la convolucion
nline = -(ntaps-1)/2:(ntaps-1)/2';
h1 = sinc(nline-delay);
ntaps=311; % Mantenerlo impar para simplificar
group_delay = (ntaps-1)/2; % Este delay hay que compensarlo porque representa la cola de la convolucion
nline = -(ntaps-1)/2:(ntaps-1)/2';
h2 = sinc(nline-delay);
tone_dig1=transpose(tone_dig);
tone_dig1=[tone_dig1; zeros(group_delay,1)];
tone_dig_delay = filter(h,1,tone_dig1);% Agrego zeros para mantener el largo de la senial filtrada
tone_dig_delay=tone_dig_delay(group_delay+1:end); % Corrijo el retardo de grupo


%% PLOT EN TIEMPO
figure
% Subplot de arriba: se muestra la se�al discreta y la se�al continua antes del delay
subplot(2, 1, 1);
p1 = plot(t_up_v/1e-6, tone_cont ,'-r','Linewidth',1);
p1.MarkerFaceColor = p1.Color;
p1.MarkerEdgeColor = p1.Color;
hold on;
p2 = plot(t_v/1e-6, tone_dig ,'ob','Linewidth',1.5);
xlabel('Tiempo (us)');
ylabel('Amplitud');
title('Se�ales discretas y continuas antes del delay');
legend([p1, p2], 'Se�al continua', 'Se�al discreta');
xlim([T_tone/1e-6, 3*T_tone/1e-6]); 
% Subplot de abajo: se muestra la se�al discreta y la se�al continua despu�s del delay
subplot(2, 1, 2);
xlabel('Tiempo (us)');
ylabel('Amplitud');
p3 = plot(t_up_v/1e-6, delay_analogico ,'--k','Linewidth',1);
p3.MarkerFaceColor = p3.Color;
p3.MarkerEdgeColor = p3.Color;
hold on;
p4 = plot(t_v/1e-6, tone_dig_delay ,'ob','Linewidth',1.5);
title('Se�ales discretas despu�s del delay ');
legend([p3, p4], 'Se�al continua despu�s del delay', 'Se�al discreta despu�s del delay');
% Mostrar solo dos periodos 
xlim([T_tone/1e-6, 3*T_tone/1e-6]);
%%


%FFT variables
NFFT = 16 * n_samples;
f_v = (-NFFT/2:NFFT/2-1)*fs/NFFT;
f_up_v = (-NFFT*OVR/2:NFFT*OVR/2-1)*fs_ch/(NFFT*OVR);

% FFTs
X_d = fftshift(abs(fft(tone_dig, NFFT))/n_samples);
X_c = fftshift(abs(fft(tone_cont, NFFT*OVR))/(n_samples*OVR));
X_d_delay = fftshift(abs(fft(tone_dig_delay, NFFT))/(n_samples));
X_c_delay = fftshift(abs(fft(delay_analogico, NFFT*OVR))/(n_samples*OVR));
H_v = fftshift(abs(fft(h, NFFT)));

% Freq
figure
plot(f_up_v/1e6, X_c,'-r','Linewidth',2);
hold on;
plot(f_up_v/1e6, X_c_delay ,'--k','Linewidth',2);
plot(f_v/1e6, X_d ,'--b','Linewidth',2);
plot(f_v/1e6, X_d_delay ,'-.c','Linewidth',1);
xlim([-201, 201]);
ylim([0, 1.3]);
plot(f_v/1e6, H_v ,'-.m','Linewidth',1);
tit = sprintf('Respuesta en frecuencia');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Frequency [MHz]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'C','C-delay','D-delay','D-ideal', 'H-Fsinc'}, 'Interpreter','latex','FontSize', fz-2);
%grid on;
% set(gcf, 'Position', [550 50 800 500],'Color', 'w');
%Plot en freq para ver si el filtro afecta el modulo de la señal
% figure
% Ls=length(tone_dig);
% TONE_DIG=fft(1/Ls.*tone_dig.*hamming(Ls), NFFT);
% fv=0:pi/NFFT:2*pi-pi/NFFT;
% sl=1:NFFT/2;
% H = fft(h,NFFT); % Filtro
% Ls=length(yf);
% YF=fft(1/Ls.*yf.*hamming(Ls), NFFT);
% 
% plot(fv(sl),abs(TONE_DIG(sl)))
% hold all
% plot(fv(sl),abs(H(sl)))
% plot(fv(sl),abs(YF(sl)))
