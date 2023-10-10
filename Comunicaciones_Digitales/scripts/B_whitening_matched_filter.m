% Programmer(s): Santiago F. Leguizamon (stgoleguizamon@gmail.com)

%                       SIM BLOCK DIAGRAM
%                       *****************
%           __________     __________     __________    
%          |          |   |          |   |          |   
%     ak---|   M(Z)   |---| M*(1/z*) |---|1/M*(1/z*)|---
%          |__________|   |__________|   |__________|   
%                          __________     __________    
%                         |          |   |          |   
%          white noise ---| M*(1/z*) |---|1/M*(1/z*)|--- white noise
%                         |__________|   |__________|   

clear all; clc; close all;

fs                  = 32e9  ;
baudrate            = 32e9  ;

n_symbs             = 1e5   ;
mod_type            = 'qam' ;
mapper_levels       = 4     ;

noise_var = 0.5;

f0 = baudrate;
w0 = pi/(fs/baudrate); % Notch frequency

% A notch filter frequency response H(z) is given by 
% H(z)=( 1 - z^-1 * e^(jw0) )*( 1 - z^-1 * e^(-jw0) )
% Impulse response is h[n] = {1, -2*cos(w0), 1}

% Let's define  M(z)=( 1 - z^-1 * e^(jw0) ) where m[n] = {1, -e^(jw0)},
%               M*(1/z*)=( 1 - z^-1 * e^(-jw0) ) where m*[n] = {1, -e^(-jw0)} and
%               gamma^2 = 1
% Thus H(z) = gamma^2.M(z).M*(1/z*)

% Local variables
bits_per_symbol = log2(mapper_levels);

% Filters
m_b_taps_v = [1, -exp(1j*w0)];
m_a_taps_v = [1];

m_conj_b_taps_v = [1, -exp(-1j*w0)];
m_conj_a_taps_v = [1];

inv_m_conj_b_taps_v = [1];
inv_m_conj_a_taps_v = [1, -exp(-1j*w0)];

h_b_taps_v = [1, -2*cos(w0), 1];
h_a_taps_v = [1];

% Processing
%bits_v = bits_generator(bits_per_symbol*n_symbs);
%symbols_v = mapper(bits_v, mapper_levels, mod_type).';

randi_v = randi([0 mapper_levels - 1], n_symbs, 1);
symbols_v = qammod(randi_v, mapper_levels);

symbols_up_v = upsample(symbols_v, fs/baudrate);

%noise_v = wgn_generator(noise_var/2, size(symbols_up_v)) + ...
%                    1j*wgn_generator(noise_var/2, size(symbols_up_v));
                
noise_v = sqrt(noise_var/2) * (randn(size(symbols_up_v))+1j*randn(size(symbols_up_v)));

% conv(symbs, m[n])
m_out_v = filter(m_b_taps_v, m_a_taps_v, symbols_up_v);

% conv(m_out, m*[n]) (signal component)
m_conj_out_s_v = filter(m_conj_b_taps_v, m_conj_a_taps_v, m_out_v);

% conv(m_out, m*[n]) (noise component)
m_conj_out_n_v = filter(m_conj_b_taps_v, m_conj_a_taps_v, noise_v);

% conv(m_conj_out, 1/m*[n]) (signal component)
inv_m_conj_out_s_v = filter(inv_m_conj_b_taps_v, inv_m_conj_a_taps_v, m_conj_out_s_v);

% conv(m_conj_out, 1/m*[n]) (noise component)
inv_m_conj_out_n_v = filter(inv_m_conj_b_taps_v, inv_m_conj_a_taps_v, m_conj_out_n_v);

%% Filters magnitude responses
figure(1); clf
[H, W] = freqz(m_b_taps_v, m_a_taps_v);
plot(W/pi*fs/2/1e9, 10*log10(abs(H)), 'DisplayName', '$M(z)$', 'LineWidth', 1.5);
hold on

[H, W] = freqz(m_conj_b_taps_v, m_conj_a_taps_v);
plot(W/pi*fs/2/1e9, 10*log10(abs(H)), '--','DisplayName', '$M^*(1/z^*)$', 'LineWidth', 1.5);

[H, W] = freqz(inv_m_conj_b_taps_v, inv_m_conj_a_taps_v);
plot(W/pi*fs/2/1e9, 10*log10(abs(H)), 'DisplayName', '$1/M^*(1/z^*)$', 'LineWidth', 1.5);
grid on

[H, W] = freqz(h_b_taps_v, h_a_taps_v);
plot(W/pi*fs/2/1e9, 10*log10(abs(H)), 'DisplayName', '$H(z)$', 'LineWidth', 1.5);

grid on
title('Filters magnitude responses', 'Interpreter', 'latex')
legend({'$M(z)$','$M^*(1/z^*)$','$1/M^*(1/z^*)$', '$Sh(z)$'},'Interpreter', 'latex')
xlabel('Frequency [GHz]', 'Interpreter', 'latex')
set(gcf, 'color', 'w', 'Position', [0 0 500 500])

%% Signal PSD
NFFT = 1024*8;
WELCH_OVERLAP = 0.*NFFT;

figure(2); clf;
[Ps, f_Ps] = pwelch(symbols_v, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs) ;
Ps_dB = 10*log10(Ps);
Ps_dB = Ps_dB - max(Ps_dB);
plot(f_Ps/1e9, Ps_dB, '-b', 'DisplayName', 'Symbols PSD')
hold on;

[Ps, f_Ps] = pwelch(m_out_v, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs) ;
Ps_dB = 10*log10(Ps);
Ps_dB = Ps_dB - max(Ps_dB);
plot(f_Ps/1e9, Ps_dB, '-g', 'DisplayName', '$M(z)$  output PSD')

[Ps, f_Ps] = pwelch(m_conj_out_s_v, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs) ;
Ps_dB = 10*log10(Ps);
Ps_dB = Ps_dB - max(Ps_dB);
plot(f_Ps/1e9, Ps_dB, '-r', 'DisplayName', '$M^*(1/z^*)$ output PSD')


[Ps, f_Ps] = pwelch(inv_m_conj_out_s_v, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs) ;
Ps_dB = 10*log10(Ps);
Ps_dB = Ps_dB - max(Ps_dB);
plot(f_Ps/1e9, Ps_dB, '--k', 'DisplayName', '$1/M^*(1/z^*)$ output PSD')
hold on;
grid on;
ylim([-70, 0.5])
xlim([min(f_Ps/1e9), max(f_Ps/1e9)])
title('Signal PSD','Interpreter', 'latex')
legend({'Symbols PSD','$M(z)$  output PSD','$M^*(1/z^*)$ output PSD','$1/M^*(1/z^*)$ output PSD'},...
                            'Interpreter', 'latex', 'Location', 'southeast')
xlabel('Frequency [GHz]', 'Interpreter', 'latex')
set(gcf, 'color', 'w', 'Position', [0 0 500 500])

%% Noise PSD
NFFT = 1024*8;
WELCH_OVERLAP = 0.*NFFT;

figure(3); clf;
[Ps, f_Ps] = pwelch(noise_v, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs) ;
Ps_dB = 10*log10(Ps);
Ps_dB = Ps_dB - max(Ps_dB);
plot(f_Ps/1e9, Ps_dB, '-b', 'DisplayName', '$N_0$')
hold on;

[Ps, f_Ps] = pwelch(m_conj_out_n_v, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs) ;
Ps_dB = 10*log10(Ps);
Ps_dB = Ps_dB - max(Ps_dB);
plot(f_Ps/1e9, Ps_dB, '--r', 'DisplayName', '$N_0|M^*(1/z^*)|^2$')

[Ps, f_Ps] = pwelch(inv_m_conj_out_n_v, hanning(NFFT/2), WELCH_OVERLAP, NFFT, fs) ;
Ps_dB = 10*log10(Ps);
Ps_dB = Ps_dB - max(Ps_dB);
plot(f_Ps/1e9, Ps_dB, '--k', 'DisplayName', '$N_0|1/M^*(1/z^*)|^2$')
hold on;
grid on;
ylim([-70, 0.5])
xlim([min(f_Ps/1e9), max(f_Ps/1e9)])
title('Noise PSD','Interpreter', 'latex')
legend({'$N_0$', '$N_0|M^*(1/z^*)|^2$', '$N_0|1/M^*(1/z^*)|^2$'}, ...
                                'Interpreter', 'latex', 'Location', 'southeast')
xlabel('Frequency [GHz]', 'Interpreter', 'latex')
set(gcf, 'color', 'w', 'Position', [0 0 500 500])
