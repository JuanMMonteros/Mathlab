clear all
close all

%% 
% Se requiere generar una señal "en tiempo continuo" formada por varios
% tonos (cos() o sin()). Luego se requiere muestrear dicha senial con un
% ADC de "resolucion infinita" (i.e., permite tomar las muestras sin
% agregar ruido de cuantizacion). Observar el espectro (usando FFT) a la
% salida del conversor para distinstas frecuencias de tono

fs = 150e6; %[Hz] frecuencia de muestreo del conversor
NOS = 16; % Oversampling factor para emular "tiempo continuo"
fch=NOS*fs; % Frecuencia del dominio de tiempo continuo
Tch=1/fch;
Ts=1/fs;

ejemplo=1;

switch ejemplo
        
        % Ejemplo 1
        f0 = 2.5e6; % Frecuencia del tono
        Tend = 100/f0; % Tiempo total de simulacion
        Nend = Tend/Tch;
        tline = (0:Nend-1).*Tch;
        tone0 = cos(2*pi*f0.*tline);
        x_t = tone0;
    case 2
        % Ejemplo 2
        f0 = 15e6; % Frecuencia del tono
        f1 = 30e6; % Frecuencia del tono
        f2 = 100e6; % Frecuencia del tono
        Tend = 100/f0; % Tiempo total de simulacion 
        Nend = Tend/Tch;
        tline = (0:Nend-1).*Tch;
        tone0 =   cos(2*pi*f0.*tline);
        tone1 = 2*cos(2*pi*f1.*tline);
        tone2 = 1.5*cos(2*pi*f2.*tline);
        x_t = tone0+tone1+tone2;
end


%% ADC
x_adc = x_t(1:NOS:end); % Muestreo la senial en TC a la fcrecuencia fs
t_adc = (0:length(x_adc)-1)*Ts; % La base de tiempo es kTs

figure
plot(tline, x_t, '-'); hold all
plot(t_adc, x_adc, 'o');
grid on
xlabel('Time [s]')
ylabel('Amplitude [V]')
legend("ADC Input", "ADC Output")

%%

% Ventaneo para evitar "spectral leakage"
W0 = hamming(length(x_t))';
W1 = hamming(length(x_adc))';
Nx = length(x_t);
NFFT = 256*1024; % Tomo una FFT muy grande para emular una TFTD (i.e., omega continua)
spectrum_tc =         1/Nx * abs(fft(W0.*x_t, NFFT)); % 1/Nx es un factor de escala que sale de calcular la FFT de un tono de duracion finita
spectrum_td = NOS  *  1/Nx * abs(fft(W1.*x_adc, NFFT)); % NOTAR EL FACTOR DE ESCALA NOS (Esto compensa el 1/T en las ecuaciones) 
fvec_tc = (0:NFFT-1)/NFFT*fch; % El vector de frecuencia llega hasta la frecuencia de muestreo de la senial en cuestion
fvec_td = (0:NFFT-1)/NFFT*fs; % El vector de frecuencia llega hasta la frecuencia de muestreo de la senial en cuestion

figure
plot(fvec_tc(1:NFFT/2), spectrum_tc(1:NFFT/2),'-'); % Dibujo las frecuencias positivas solamente
hold all
plot(fvec_td(1:NFFT/2), spectrum_td(1:NFFT/2),'--','LineWidth',2);
grid on
xlabel('Frequency [Hz]')
ylabel('Amplitude [V]')
legend("Spectrum ADC Input", "Spectrum ADC Output")
xlim([0,fs])