clear all
close all

rand('seed',1);
randn('seed',51);

% Basic TX QPSK o QAM16
BR = 32e9; %Bd
L = 10000; % Simulati on Length
fs = BR; % Sampling rate to emulate analog domain
T = 1/BR; % Time interval between two consecutive symbols
Ts = 1/fs; % Time between 2 conseutive samples at Tx output
M=16;
EbNo_dB = 15; % dB

% Parametros de la portadora
delta_freq = 10e6; % Offset del LO
LW = 00e3; % Ancho de linea [Hz]
theta0 = 0/180*pi; % Fase inicial
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
noise_real = sqrt(No/2)*randn(size(xsymb)); % Lo que acompaña al ruido es el sigma (desvio estandard)
noise_imag = sqrt(No/2)*randn(size(xsymb));
noise = noise_real + 1j*noise_imag;
% NO VALE noise= (1+1j)*randn(size(yup))
rx_noisy = xsymb + noise;
clear noise noise_real noise_imag

%
% Agregar efectos de portadora
Ldata= length(rx_noisy);
time = (0:Ldata-1).'.*Ts;

lo_offset = exp(1j*2*pi*delta_freq*time); % LO OFFSET
phase_offset = exp(1j*theta0); % STATIC PHASE
% Fluctuaciones de frecuencia
freq_fluctuations = exp(1j*frequency_fluctuations_amp/frequency_fluctuations_freq.*sin(2*pi*frequency_fluctuations_freq.*time));
% Para el ruidod de fase
freq_noise = sqrt(2*pi*LW/fs).*randn(Ldata,1);
phase_noise = cumsum(freq_noise); % Proceso de Wiener
osc_pn = exp(1j.*phase_noise);
phase_tone = exp(1j.*phase_tone_amp.*sin(2*pi*phase_tone_freq.*time));
rxs = rx_noisy.*lo_offset.*phase_offset.*osc_pn.*freq_fluctuations.*phase_tone;

scatterplot(rxs)

% figure
% plot(phase_noise)
% grid on
% xlabel('Time')
% ylabel('Channel Phase Noise [rad]')

%% PLL RX
Ldata = length(rxs);
Kp = 15e-2;
Ki = Kp/1000;
phase_error = zeros(Ldata,1);
nco_output = zeros(Ldata,1);
integral_branch = zeros(Ldata,1);
pll_output = zeros(Ldata,1);

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

figure
plot(phase_error)
grid on
xlabel("Time")
ylabel("Phase Error")

figure
plot(integral_branch*BR/2/pi/1e6)
grid on
xlabel("Time")
ylabel("Integral Branch [MHz]")

figure
plot(nco_output); hold all
plot(phase_noise); hold all
grid on
xlabel("Time")
ylabel("Phase [rad]")
legend("PLL NCO","Channel")

scatterplot(pll_output(Ldata/2:end))


%% Slicer

z = pll_output;
Ldata = length(z);
a_hat = zeros(Ldata,1);
for m=1:Ldata
    a_hat(m) = slicer(z(m), M);
end
figure
plot(real(z), '.r');
hold all
plot(real(a_hat), 'xb');
legend("Slicer Input","Slicer Output")

%% Calcular la BER
% 1. Alinear
d0 = finddelay(xsymb, z);
guard0 = fix(Ldata*.25);
guard1 = 1e3;
a_hat_align = a_hat (1+d0+guard0:end-guard1);
pll_output_align = pll_output (1+d0+guard0:end-guard1);
x_align = xsymb (1+guard0:end-guard1-d0);

% 2. OPCION 1 (millenial, facil)
simbolos_totales = length(x_align);
simbolos_errados = sum(x_align ~= a_hat_align);
symbol_error_rate = simbolos_errados / simbolos_totales;
aprox_ber_sin_corregir_cs = 1/log2(M) * symbol_error_rate

% %%
% bits_tx = demodulador_M_QAM(x_align, M);
% bits_rx = demodulador_M_QAM(a_hat_align, M);
% [errors, ber] = biterr(bits_tx, bits_rx);
% ber;


% 
% %% Correccion dinamica de los cycle slips
% Hacerlo con los simbolos alineados
% orx son los simbolos a la salida del BPS (sin slicer)
% otx son los simbolos transmitidos
% orx y otx estan alineados
orx = pll_output_align;
otx = x_align;
WINDOW_LEN = 500;
Ldata = length(orx);
nblocks = fix(Ldata/WINDOW_LEN);

orx_cs_fixed = zeros(nblocks*WINDOW_LEN,1);
cs_phase = zeros(nblocks,1);
cs_count =0;
last_phase=0;
mse_log = zeros(nblocks, 4);
phase_test_list = [0, pi/2, -pi/2, pi];

for nblock=1:nblocks
    slice = (nblock-1)*WINDOW_LEN+1: nblock*WINDOW_LEN;
    rx_block_in = orx(slice);
    tx_block_in = otx(slice);
    
    min_mse = inf;
    phase_ok = 0;
    for np = 1:4
        phase_test = phase_test_list(np);
        block_test = rx_block_in.*exp(1j*phase_test);
        mse = mean(abs(block_test-tx_block_in).^2);
        mse_log(nblock, np) = mse;
        if mse < min_mse
            min_mse = mse;
            phase_ok=phase_test;
        end
    end
    if phase_ok ~= last_phase
        cs_count = cs_count +1;
    end
    last_phase = phase_ok;
    cs_phase(nblock) = phase_ok; 
    orx_cs_fixed(slice) = rx_block_in.*exp(1j*phase_ok);
end

% figure
% plot(cs_phase/(pi/2))
% ylim([-2,2])

Ldata = length(orx_cs_fixed);
orx_cs_fixed_slicer = zeros(Ldata,1);
for m=1:Ldata
    orx_cs_fixed_slicer(m) = slicer(orx_cs_fixed(m), M);
end

simbolos_totales = length(orx_cs_fixed_slicer);
simbolos_errados = sum(x_align ~= orx_cs_fixed_slicer);
symbol_error_rate = simbolos_errados / simbolos_totales;
aprox_ber_corrigiendo_cs = 1/log2(M) * symbol_error_rate;

ber_teo = berawgn(EbNo_dB, 'qam', M)