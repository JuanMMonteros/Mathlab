clear all
close all

Lsim=4*1024;
xsig = get_limited_band_signal(Lsim,0.5*pi);

% FFT de la senial original
figure
NFFT=4*1024;
Ls=length(xsig);
XSIG=fft(1/Ls.*xsig.*hamming(length(xsig)), NFFT);
fv=0:pi/NFFT:2*pi-pi/NFFT;
sl=1:NFFT/2;
plot(fv(sl),abs(XSIG(sl)))
grid on
xlabel('Frecuencia Discreta')
ylabel('Amplitud')
title('TF de la senial original')
% Se pretende retrasar la senial una fraccion de muestra
delay=0.5; % esto representa una fraccion de muestra
% Genero el filtro de la eq 4.65 del libro
ntaps=31; % Mantenerlo impar para simplificar
group_delay = (ntaps-1)/2; % Este delay hay que compensarlo porque representa la cola de la convolucion
nline = -(ntaps-1)/2:(ntaps-1)/2;
h1 = sinc(nline-delay);

figure
stem(nline, h1);
grid on
xlabel('Discrete time')
title('Interpolation filter')

% Version alternativa del sinc
h2= my_rcosine(1,1,0.2,ntaps,delay);

yf = filter(h2,1,[xsig; zeros(group_delay,1)]); % Agrego zeros para mantener el largo de la senial filtrada
yf=yf(group_delay+1:end); % Corrijo el retardo de grupo

% Plot en tiempo discreto para ver el retardo
figure
plot(xsig,'-o');
hold all
plot(yf,'-x');

% Plot en freq para ver si el filtro afecta el modulo de la señal
figure
NFFT=4*1024;
Ls=length(xsig);
XSIG=fft(1/Ls.*xsig.*hamming(Ls), NFFT);
fv=0:pi/NFFT:2*pi-pi/NFFT;
sl=1:NFFT/2;
H1 = fft(h1,NFFT); % Filtro
H2 = fft(h2,NFFT); % Filtro
Ls=length(yf);
YF=fft(1/Ls.*yf.*hamming(Ls), NFFT);

plot(fv(sl),abs(XSIG(sl)))
hold all
plot(fv(sl),abs(H1(sl)))
plot(fv(sl),abs(H2(sl)))
%plot(fv(sl),abs(YF(sl)))




