clear all
%close all

%% Creo un filtro target
h = fir1(6,0.7); h = h/sum(h);

% figure; stem(h)
% figure; freqz(h)

%% Filtrar una entrada random con h
% la entrada podria no ser un ruido blanco
alpha = 0;
xk = randn(500e3,1); % Datos para la adaptacion (gaussiano)
% xk = filter(alpha, [1 -(1-alpha)], xk);
xk = filter(1, [1 -alpha], xk);
xk = xk/std(xk); % Normalizo para tener varianza 1


fprintf("Varianza de la señal de entrada: %2.2e", var(xk))
L =1;
corrv = xcorr(xk,xk,L, 'biased'); % Calculo correlacion considerando 1 tap de memoria
r = corrv(L+1:end); % Me quedo con la porcion causal de la auto correlacion
CorrMtx = toeplitz(r, conj(r)); % Calculo la matriz como esta en el libro

CorrMtx

eig(CorrMtx)
%%
zk = filter(h,1,xk);
nk = 0.1*randn(length(zk),1);
yk = zk + nk;
%%
NTAPS = 15;
xbuffer = zeros(NTAPS,1);
ck = zeros(NTAPS,1); % Solucion inicial (podria ser un impulso)
beta = .25e-2/4;
gradient = zeros(NTAPS,1);

yk_prima = zeros(length(yk),1);
Ek = zeros(length(yk),1);
ck_log = zeros(length(yk), NTAPS);

INPUT_DELAY=0;
OUTPUT_DELAY = 0;

for k=INPUT_DELAY+OUTPUT_DELAY+1:length(xk)
    xbuffer(2:end) = xbuffer(1:end-1);
    xbuffer(1) = xk(k-INPUT_DELAY);
    
    %yk_prima(k) = sum(xbuffer.*ck);
    yk_prima(k) = (ck.') * xbuffer;
    
    Ek(k) = yk_prima(k) - yk(k-OUTPUT_DELAY); % Salida del filtro adap - la senial deseada
    
    gradient = Ek(k).*conj(xbuffer);
    ck_log(k,:) = ck;
    ck = ck - beta*gradient;
    
end


%% Coeficientes finales
figure
Lh = length(h); Mh = (Lh-1)/2;
time_h = -Mh:Mh;
stem(time_h, h)
hold all

%ck_eff = [zeros(INPUT_DELAY,1); ck];
ck_eff  = ck;
gd=grpdelay(ck_eff); gd0 = gd(1);
time_c= -gd0:(length(ck_eff)-gd0-1);
stem(time_c,  ck_eff)
grid on

legend('Respuesta del canal', 'Respuesta del LMS')

% figure
% freqz(ck);

%% Error instantaneo
figure
plot(Ek)
title('Error instantaneo')
grid on
xlabel('Tiempo [samples]')

%% Calculo del MSE
figure
mse = abs(Ek).^2;
mse = filter(ones(1e3,1)./1e3, 1, mse);
plot(10*log10(mse))
title('MSE')
grid on
xlabel('Tiempo [samples]')
ylabel('MSE [dB]')

%% ploteo evolucion de los coeficientes
figure
plot(ck_log)
grid on
xlabel('Tiempo [samples]')
ylabel('Magnitud de los coeficientes')
