clear all
%close all
% 
% fcut=0.2;
% orden=16;
% hn=fir1(orden, fcut,'low');
% En la funcion FIR1, fcut va de 0 a <1, donde 0 representa 0 y 1
% representa pi (fdisc) = fN (fcont)
% [H,f] = freqz(hn);
% plot(f, abs(H))
% hold all

%%
% Tengo un tono en tiempo continuo de 100MHz y otro de 150MHz
% Quiero hacer un DSP que muestree los dos tonos y filtre el tono de mas
% alta frecuencia

fs = 600e6; %[Hz] frecuencia de muestreo del conversor
NOS = 16; % Oversampling factor para emular "tiempo continuo"
fch=NOS*fs; % Frecuencia del dominio de tiempo continuo
Tch=1/fch;
Ts=1/fs;

% Genero la señal en tiempo continuo
f0 = 100e6; % Frecuencia del tono
f1 = 150e6; % Frecuencia del tono
Tend = 100/f0; % Tiempo total de simulacion
Nend = Tend/Tch;
tline = (0:Nend-1).*Tch;
tone0 = 1.5*cos(2*pi*f0.*tline);
tone1 = 0.7*cos(2*pi*f1.*tline);
x_t = tone0+tone1;

% Paso por el conversor C/D
x_adc = x_t(1:NOS:end); % Muestreo la senial en TC a la fcrecuencia fs
t_adc = (0:length(x_adc)-1)*Ts; % La base de tiempo es kTs

% Filtro
fcut_analogica=(f0+f1)/2;
fcut_digital = fcut_analogica / (fs/2);
orden=128;
hn=fir1(orden, fcut_digital,'low');
y_filtro=filter(hn,1,x_adc);

% Truquitos para dibujar bien SI HUBIERA DIBUJADO LAS SE;ALES EN TIEMPO
y_filtro=y_filtro(1+orden/2:end);
x_adc=x_adc(1:end-orden/2);
t_adc=t_adc(1:end-orden/2);

% Dibujos
NFFT = 256*1024;
[Hf,f] = freqz(hn, 1,NFFT/2);

W0 = hamming(length(x_t))';
W1 = hamming(length(x_adc))';
Nx = length(x_t);
 % Tomo una FFT muy grande para emular una TFTD (i.e., omega continua)
spectrum_tc =         1/Nx * abs(fft(W0.*x_t, NFFT)); % 1/Nx es un factor de escala que sale de calcular la FFT de un tono de duracion finita
spectrum_x_td = NOS  *  1/Nx * abs(fft(W1.*x_adc, NFFT)); % NOTAR EL FACTOR DE ESCALA NOS (Esto compensa el 1/T en las ecuaciones) 
spectrum_y_td = NOS  *  1/Nx * abs(fft(W1.*y_filtro, NFFT)); 
fvec_tc = (0:NFFT-1)/NFFT*fch; % El vector de frecuencia llega hasta la frecuencia de muestreo de la senial en cuestion
fvec_td = (0:NFFT-1)/NFFT*fs; % El vector de frecuencia llega hasta la frecuencia de muestreo de la senial en cuestion

figure
plot(fvec_tc(1:NFFT/2), spectrum_tc(1:NFFT/2),'-'); % Dibujo las frecuencias positivas solamente
hold all
plot(fvec_td(1:NFFT/2), spectrum_x_td(1:NFFT/2),'--','LineWidth',2);
plot(fvec_td(1:NFFT/2), abs(Hf),'--k','LineWidth',2);
plot(fvec_td(1:NFFT/2), spectrum_y_td(1:NFFT/2),'-','LineWidth',2);

grid on
xlabel('Frequency [Hz]')
ylabel('Amplitude [V]')
legend("Spectrum ADC Input", "Spectrum ADC Output", "Freq. Resp. Filter", "Spectrum Filter Output")
xlim([0,fs])

