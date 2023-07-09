clear all
close all

%% Ejemplo
% generar un tono en un procesador digital (la frecuencia del tono esta expresada en [0,pi])
% Modelar un DAC para transformar el tono discreto en un señal de tiempo
% continuo. La frecuencia de muestreo del DAC es 200MHz. 

fs = 200e6; % DAC sampling rate
NOS = 4; % Oversampling factor para emular "tiempo continuo"

fch=NOS*fs; % Frecuencia del dominio de tiempo continuo
Tch=1/fch;
Ts=1/fs;

Nsamples = 200; % Numero de muestras para del tono
omega_0 = 0.3*pi; % Corresponde a w0=0.1*pi*fs --> f0 = 0.1*fs/2;
A0=2;
tone0_disc = A0*cos(omega_0.*(0:Nsamples-1));
%omega_1 = 0.2*pi;
%tone1_disc =   2*cos(omega_1.*(0:Nsamples-1));
%x_disc = tone0_disc+tone1_disc;
x_disc = tone0_disc;

 %Genero una señal de referencia que es el tono que "en tiempo continuo"
tline_continuo=(0:Nsamples*NOS-1).*Tch;
tone_ref = A0*cos(omega_0*fs.*tline_continuo');

% Modelo del DAC
% Un DAC consiste en "retener" un voltage constante para cada muestra que
% ingresa en el dominio discreto. Esto se puede hacer un filtro llamado
% "retenedor de orden cero" o Zero Order Hold (ZOH). Es un filtro cuya
% respuesta al impulso es un pulso de duracion Ts

% Opcion 1: retenedor de orden 0
%Si Ts = NOS*Tch, entonces el filtro es un pulso rectangular de NOS
% muestras
dac_filter = ones(1,NOS); grpd=0;

% Opcion 2: retenedor de orden 1 (interp lineal)
%dac_filter=1/NOS*conv(ones(1,NOS),ones(1,NOS)); grpd=NOS-1;

% Opcion 3: Reconstruccion ideal
%nfilt=-1000:1000; dac_filter = sinc(nfilt/Ts*Tch); grpd = grpdelay(dac_filter);grpd =grpd(1);

% opcion 4: Filtro FIR custom
%Nord=16; dac_filter=NOS*fir1(Nord,fs/fch); grpd=Nord/2;



% TAREA: calcular la respuesta en frecuencia de un ZOH

% El ZOH esta funcionando a una tasa NOS veces la señal discreta
% Entonces, lo primero que hay que hacer es "sobremuestrear" la señal
% discreta para ajustar su tasa de muestreo a la misma tasa que el filtro
% retenedor
x_up = upsample(x_disc, NOS); % Agrega NOS-1 ceros entre cada muestra

% Finalmente, filtro
x_tc1 = filter(dac_filter,1,x_up);
x_tc = [x_tc1(grpd+1:end) zeros(1, grpd)];

% Plot en tiempo
figure
tline_disc = (0:Nsamples-1).*Ts;
tline_cont = (0:Nsamples*NOS-1).*Tch;
plot(tline_disc, x_disc,'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k');
hold all
stem(tline_cont, x_up,'x');
plot(tline_cont, x_tc,'-', 'LineWidth',2);
plot(tline_cont, tone_ref,'--k', 'LineWidth',1.5);

xlabel('Time [s]')
ylabel('Amplitude [V]')
grid on
legend("Discrete samples", "ZOH input","ZOH output")

%%
% Plot en frecuencia
NFFT = 256*1024;%resolucion de la w

fvec_tc = (0:NFFT-1)/NFFT*fch;
fvec_td = (0:NFFT-1)/NFFT*fs;

W0 = hamming(length(x_disc))';
W1 = hamming(length(x_up))';
Ns = length(x_disc);
spect_disc = 1/Ns*abs(fft(W0.*x_disc, NFFT));
spect_x_up = 1/Ns*abs(fft(W1.*x_up, NFFT));
spect_tc =   1/Ns*1/NOS*abs(fft(W1.*x_tc, NFFT));
dac_response = abs(fft(dac_filter,NFFT));

figure
subplot 311
plot(fvec_td(1:NFFT/2), spect_disc(1:NFFT/2),'-b','LineWidth',2.5);
hold all
grid on
xlabel('Frequency [Hz]');
ylabel('Amplitude [V]');
legend('Discrete sequence spectrum');

subplot 312
plot(fvec_tc(1:NFFT/2), spect_x_up(1:NFFT/2),'-b','LineWidth',3);
hold all
plot(fvec_tc(1:NFFT/2), dac_response(1:NFFT/2)/dac_response(1),'--k','LineWidth',2);
grid on
xlabel('Frequency [Hz]');
ylabel('Amplitude [V]');
legend('Upsampled discrete sequence', 'DAC Filter response');

subplot 313
plot(fvec_tc(1:NFFT/2), spect_tc(1:NFFT/2),'-b','LineWidth',3);
grid on
xlabel('Frequency [Hz]');
ylabel('Amplitude [V]');
legend('Cont. time signal' );


