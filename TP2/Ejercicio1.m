% Parámetros de configuración
fs = 200e6; % Frecuencia de muestreo del DAC
Ts=1/fs; %Periodo de Mustreo
OVR = 10; % Oversampling factor
fch = OVR * fs; % Frecuencia de muestreo del tiempo continuo
Tch=1/fch; %Periodo del canal 

Nsamples = 201; % Número de muestras para la señal discreta
omega_0 = 0.1 * pi; % Frecuencia del tono discreto
T_tone=2*pi/(omega_0*fs);% Periodo del Tono 
A0 = 2; % Amplitud del tono discreto
% Generar señal discreta cosenoidal D/C
t_disc = (0:Nsamples-1).*Ts;
x_disc = A0 * cos(omega_0.*t_disc*fs);

% Generar señal de tiempo continuo cosenoidal de referencia
t_cont = (0:Nsamples*OVR-1).*Tch;
tone_ref = A0 * cos(omega_0 * fs * t_cont);


filtro=1;

switch filtro
    case 1
        % Retenedor de orden cero (ZOH)
        dac_filter = ones(1,OVR);
        grpd=0;
    case 2
        % Interpolación lineal
       dac_filter=1/OVR*conv(ones(1,OVR),ones(1,OVR)); 
       grpd=OVR-1;
    case 3
        % Reconstrucción ideal (filtro sinc)
        nfilt=-1000:1000;
        dac_filter = sinc(nfilt/Ts*Tch); 
        grpd = grpdelay(dac_filter);
        grpd =grpd(1);
    case 4
        % Filtro FIR personalizado
        Nord = 16;
        dac_filter =OVR * fir1(Nord, fs / fch);
end

x_up = upsample(x_disc, OVR); % Agrega OVR-1 ceros entre cada muestra

% filtro
x_tc1 = filter(dac_filter,1,x_up);
x_tc = [x_tc1(grpd+1:end) zeros(1, grpd)];

%%
%Plote del tiempo 
figure
hold all
plot(t_disc*1e6, x_disc,'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k');
plot(t_cont*1e6, x_tc,'-', 'LineWidth',2);
plot(t_cont*1e6, tone_ref,'--k', 'LineWidth',1.5);
xlim([0 2*T_tone*1e6])
xlabel('Time [$\mu$s]', 'Interpreter','latex','FontSize', 12);
ylabel('Amplitude [V]', 'Interpreter','latex','FontSize', 12);
grid on
legend("Discrete samples", "Salida D/C","signal Reference");

%%
% Plot en frecuencia
NFFT = 256*1024;%resolucion de la w
fvec_tc = (0:NFFT-1)/NFFT*fch;
fvec_td = (0:NFFT-1)/NFFT*fs;

% Ventaneo para evitar "spectral leakage"
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


