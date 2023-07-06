clear all
close all

%% 
% Se requiere generar una señal "en tiempo continuo" formada por varios
% tonos (cos() o sin()). Luego se requiere muestrear dicha senial con un
% ADC de "resolucion infinita" (i.e., permite tomar las muestras sin
% agregar ruido de cuantizacion). Observar el espectro (usando FFT) a la
% salida del conversor para distinstas frecuencias de tono

fs = 2.5*200e6; %[Hz] frecuencia de muestreo del conversor
NOS = 16; % Oversampling factor para emular "tiempo continuo"
fch=NOS*fs; % Frecuencia del dominio de tiempo continuo
Tch=1/fch;
Ts=1/fs;

f0 = 150e6; % Frecuencia del tono
f1 = 30e6; % Frecuencia del tono
f2 = 78e6; % Frecuencia del tono
Tend = 100/f0; % Tiempo total de simulacion 
Nend = Tend/Tch;
tline = (0:Nend-1).*Tch;
tone0 =   cos(2*pi*f0.*tline);
tone1 = 2*cos(2*pi*f1.*tline);
tone2 = 1.5*cos(2*pi*f2.*tline);

z_t = tone0;%+tone1+tone2; % senial TC

h_aaf = fir1(64, (fs/2)/(fch/2)); % Dejo los "/2" para mayor claridad. La frecuencia de corte debe ser fN=fs/2.
x_t = filter(h_aaf,1,z_t);


%% ADC
x_adc = x_t(1:NOS:end); % Muestreo la senial en TC a la fcrecuencia fs
t_adc = (0:length(x_adc)-1)*Ts; % La base de tiempo es kTs


% Ventaneo para evitar "spectral leakage"
W0 = hamming(length(x_t))';
W1 = hamming(length(x_adc))';
Nx = length(x_t);
NFFT = 256*1024; % Tomo una FFT muy grande para emular una TFTD (i.e., omega continua)
spectrum_tc_before_aaf = 1/Nx * abs(fft(W0.*z_t, NFFT)); % 1/Nx es un factor de escala que sale de calcular la FFT de un tono de duracion finita
spectrum_tc_after_aaf =  1/Nx * abs(fft(W0.*x_t, NFFT)); % 1/Nx es un factor de escala que sale de calcular la FFT de un tono de duracion finita
spectrum_td =    NOS  *  1/Nx * abs(fft(W1.*x_adc, NFFT)); % NOTAR EL FACTOR DE ESCALA NOS (Esto compensa el 1/T en las ecuaciones) 
HAAF = abs(fft(h_aaf, NFFT));

fvec_tc = (0:NFFT-1)/NFFT*fch; % El vector de frecuencia llega hasta la frecuencia de muestreo de la senial en cuestion
fvec_td = (0:NFFT-1)/NFFT*fs; % El vector de frecuencia llega hasta la frecuencia de muestreo de la senial en cuestion

figure
plot(fvec_tc(1:NFFT/2), spectrum_tc_before_aaf(1:NFFT/2),'-b', 'LineWidth',1.5); % Dibujo las frecuencias positivas solamente
hold all
plot(fvec_tc(1:NFFT/2), spectrum_tc_after_aaf(1:NFFT/2),'-g', 'LineWidth',2); % Dibujo las frecuencias positivas solamente
plot(fvec_td(1:NFFT/2), spectrum_td(1:NFFT/2),'--r','LineWidth',2);
plot(fvec_tc(1:NFFT/2),HAAF(1:NFFT/2),'--k', 'LineWidth',2);
grid on
xlabel('Frequency [Hz]')
ylabel('Amplitude [V]')
legend("Spectrum AAF Input", "Spectrum ADC Input", "Spectrum ADC Output","AAF Response")
xlim([0,fs/2+50e6])
