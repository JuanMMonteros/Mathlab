%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero
%-----------------------------------------------------------------------------%

clear; close all; clc;

%% Basic TX QPSK 

M = 4;          % QPSK
L = 1e5;        % Simulation Length
BR = 32e9;      % Baud Rate
N = 4;          % Oversampling rate
rolloff = 0.5;  % Pulse shaping rolloff
h_taps = 101;   % Pulse shaping taps
EbNo_db = 5;    % EbNo in dB

fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 consecutive samples at Tx output

% Two symbols generation (+1,-1) for QPSK
ai = 2*randi([0,1],L,1)-1; % Real
aq = 2*randi([0,1],L,1)-1; % Imag
ak = ai + 1j*aq;

% Upsampling to change sampling rate
xup = upsample(ak,N);

% Filter to interpolate the signal
h = root_raised_cosine(BR/2, fs, rolloff, h_taps, 0);
h_delay = (h_taps-1)/2;
yup = filter(h,1,[xup; zeros(h_delay, 1)]);
yup = yup(1+h_delay:end);

%% Channel

% EbNo to channel snr
k = log2(M);
EbNo = 10^(EbNo_db/10);
SNR_slc = EbNo * k;
SNR_ch = SNR_slc / N;

% Noise power
Ps = var(yup);
Pn = Ps/SNR_ch;

% Noise generator
n = sqrt(Pn/2) .* (randn(length(yup),1) + 1j.*randn(length(yup),1));

% Noise addition
rx = yup + n;

%% RX

% Filter with MF
h_mf = conj(h(end:-1:1));
ymf = filter(h_mf,1,[rx; zeros(h_delay, 1)]);
ymf = ymf(1+h_delay:end);

% Downsampling
rx_down = ymf(1:N:end);

% Constellation Plot
scatterplot(rx_down(100:10000));
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

%% Slicer and BER estimator

% Slicer
ak_hat = my_slicer(rx_down, M);

% Theo ber
ber_theo = berawgn(EbNo_db, 'qam', M);

% Estimated ber
[ber_sim, n_errors] = my_ber_checker(ak_hat, ak, M, 'auto');

% Print
fprintf('\t\n - Theo BER = %.2e\n', ber_theo)
fprintf('\t\n - Sim BER = %.2e\n', ber_sim)
fprintf('\t\n - Errors = %d\n', n_errors)