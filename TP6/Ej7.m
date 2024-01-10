clear all
close all

rand('seed',1);
randn('seed',51);

% Basic TX QPSK o QAM16
position_and_size = [100 100 800 600];
BR = 32e9; %Bd
L = 10000; % Simulati on Length
fs = BR; % Sampling rate to emulate analog domain
T = 1/BR; % Time interval between two consecutive symbols
Ts = 1/fs; % Time between 2 conseutive samples at Tx output
M=16;
EbNo_dB = 10; % dB

% Parametros de la portadora
delta_freq = 10e6; % Offset del LO
LW = 00e3; % Ancho de linea [Hz]
theta0 = 15/180*pi; % Fase inicial
frequency_fluctuations_amp = 0e6;
frequency_fluctuations_freq = 1e3;
phase_tone_amp = 10/180*pi;
phase_tone_freq = 250e6;

% Calculos de la SNR
EbNo = 10^(EbNo_dB/10);
SNR_Slc = EbNo * log2(M);
SNR_canal = SNR_Slc;

% Two symbols generation (+1,-1) for QPSK
xsymb = qammod(randi([0 M-1], L,1), M);

% Agrego ruido
% SNR_canal [dB] = EbNo [dB] + 10*log10(k) - 10*log10(N); % Formula que por ahi se ve en internet
% k = log2(M)
Psignal = var(xsymb);
Pnoise = Psignal/SNR_canal;
No = Pnoise; % varianza por unidad de frecuencia
noise_real = sqrt(No/2)*randn(size(xsymb)); % Lo que acompaÃ±a al ruido es el sigma (desvio estandard)
noise_imag = sqrt(No/2)*randn(size(xsymb));
noise = noise_real + 1j*noise_imag;
% NO VALE noise= (1+1j)*randn(size(yup))
rx_noisy = xsymb + noise;
clear noise noise_real noise_imag

%
%% PLL RX
Ldata = length(rx_noisy);
Kp = 0.05;
Ki = Kp/1000;
phase_error = zeros(Ldata,1);
nco_output = zeros(Ldata,1);
integral_branch = zeros(Ldata,1);
pll_output = zeros(Ldata,1);
% Agregar efectos de portadora
time = (0:Ldata-1).'.*Ts;
frec =logspace(log10(1e0), log10(fs/2), 1024); % Barrido de frec.
pll_out = zeros(length(frec),1);
for i = 1:length(frec)
    % -- Barrido de senoidales -- %
    fn = frec(i);
    wn = 2*pi*fn;
    tita0 = 15;
    tita_in = tita0/180*pi * sin(wn.*time);

rxs = rx_noisy.*exp(1j.*tita_in);



        for m=2:Ldata
            xpll = rxs(m);
            derot_x = xpll*exp(-1j*nco_output(m-1));
            pll_output(m) = derot_x;
            a_hat = slicer(derot_x, M);
            phase_error(m) = angle(derot_x .* conj(a_hat));
            prop_error = phase_error(m) * Kp;
            integral_branch(m) = integral_branch(m-1) + Ki*phase_error(m);
            nco_output(m) = nco_output(m-1) + prop_error +integral_branch(m);
        end
        pll_out(i) = 10*log10(var(nco_output)/var(tita_in));
end

n_freq_pos = 2^16;  
n_freq_frac = 2^16;
n_v = [0, logspace(-log10(n_freq_frac),log10(n_freq_pos),n_freq_pos+n_freq_frac-1)];
f_v = n_v * BR/2/n_freq_pos;
wd_v = n_v * pi/n_freq_pos;

G_th_v = (Kp+Ki.*1./(1-exp(-1j*wd_v))) * 1./(1-exp(-1j.*wd_v));
H_th_v = G_th_v ./ (1+G_th_v.*exp(-1j*wd_v));
H_th_db_v = 20*log10(abs(H_th_v));
 % Eliminar la última muestra de pll_output y frec
pll_out = pll_out(1:end-1);
frec = frec(1:end-1);


figure
semilogx(frec, pll_out,'-r', 'LineWidth', 2);
hold all 
semilogx(f_v, H_th_db_v, '--b', 'LineWidth', 2);
grid on 
xlabel('Frecuencia [Hz]');
ylabel('Amplitude [dB]');
xlim([2e6 4e8]);
title('Jitter Transfer Plot');
xlabel('Time [s]')
ylabel('Amplitudes [dB]')
legend('Simulado', 'Teórico'); % Agrega la leyenda
set(gcf, 'Position',position_and_size,'Color', 'w');


%% Slicer

z = pll_output;
Ldata = length(z);
a_hat = zeros(Ldata,1);
for m=1:Ldata
    a_hat(m) = slicer(z(m), M);
end
