%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% The nomenclature "_v" at the end of a variable indicates that is a vector
% The nomenclature "_m" at the end of a variable indicates that is a matrix
% The nomenclature "_s" at the end of a variable indicates that is a struct
% The nomenclature "_c" at the end of a variable indicates that is a cell

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

%% Parameters

fs = 1e6;       % Sampling frequency
ts = 1/fs;      % Sampling time

OVR = 10;       % Tone oversampling
f0 = fs/OVR;    % Tone frequency
T0 = 1/f0;      % Tone period

P = 50;          % Number of periods for the tone

%% Processing

% Time vector
n_samples = P*OVR;
t_v = (0:n_samples-1)*ts;

% -- Signals creation --
x_v = sin(2*pi*f0*t_v);
n_v = 1*randn(1, n_samples);
x_noisy_v = x_v + n_v;

% IIR filter creation
order = 3;
fc1 = 0.9 * f0 / (fs/2);
fc2 = 1.1 * f0 / (fs/2);
[b_iir, a_iir] = butter(order, [fc1 fc2], 'bandpass');

% I create a FIR that aproximates the IIR to plot the response
impulse_taps = 200;
fir_delay = impulse_taps/2;
impulse_v = [1, zeros(1,impulse_taps-1)];
h_fir_v = filter(b_iir, a_iir, impulse_v);
%figure; plot(h_fir_v); % Here I see that the impulse length is enough 

% Filter the signal
x_fil_v = filter(b_iir,a_iir,x_noisy_v);

% -- Spectrum --
NFFT = 10*n_samples;
X_v = fft(x_v,NFFT)./n_samples;
X_noisy_v = fft(x_noisy_v,NFFT)./n_samples;
X_fil_v = fft(x_fil_v,NFFT)./n_samples;
H_v = fft(h_fir_v,NFFT);

% Frequency vector
f_v = (-NFFT/2:NFFT/2-1)*fs/NFFT;

%% Plots
fz = 15;

figure;

% Signals in time
subplot(2,2,1);
plot(t_v/1e-6,x_v,'-r','Linewidth',2)
hold on;
plot(t_v/1e-6,x_noisy_v,'--b','Linewidth',2)
plot(t_v/1e-6,x_fil_v,'-.m','Linewidth',2)

title('Signals in time','Interpreter','latex','FontSize', fz);
leg_c = {'Ideal','Noisy','Filtered'};
legend(leg_c,'Location','ne','Interpreter','latex','FontSize', fz-2);
xlabel('Time [useg]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
grid on

% h(n)
subplot(2,2,3);
stem(h_fir_v,'r','Linewidth',2)
hold on;
title('h(n)','Interpreter','latex','FontSize', fz);
xlabel('Samples', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
grid on

% Module in freq
subplot(2,2,2);
plot(f_v/1e3,fftshift(abs(X_v)),'-r','Linewidth',2)
hold on;
plot(f_v/1e3,fftshift(abs(X_noisy_v)),'--b','Linewidth',2)
plot(f_v/1e3,fftshift(abs(X_fil_v)),'-.m','Linewidth',2)
plot(f_v/1e3,fftshift(abs(H_v)),'-k','Linewidth',1)
title('Module in frequency','Interpreter','latex','FontSize', fz);
leg_c = {'Ideal','Noisy','Filtered','abs(H)'};
legend(leg_c,'Location','nw','Interpreter','latex','FontSize', fz-2);
xlabel('Frequency [KHz]', 'Interpreter','latex','FontSize', fz);
ylabel('Module', 'Interpreter','latex','FontSize', fz);
grid on; 

% Phase in freq
subplot(2,2,4);
plot(f_v/1e3,fftshift((angle(X_v))),'-r','Linewidth',2)
hold on;
plot(f_v/1e3,fftshift(angle(X_noisy_v)),'--b','Linewidth',2)
plot(f_v/1e3,fftshift(angle(X_fil_v)),'-.m','Linewidth',2)
plot(f_v/1e3,fftshift(angle(H_v)),'-k','Linewidth',1)
title('Angle in frequency','Interpreter','latex','FontSize', fz);
leg_c = {'Ideal','Noisy','Filtered','angle(H)'};
legend(leg_c,'Location','nw','Interpreter','latex','FontSize', fz-2);
xlabel('Frequency [KHz]', 'Interpreter','latex','FontSize', fz);
ylabel('Angle', 'Interpreter','latex','FontSize', fz);
grid on; 

set(gcf, 'Position', [50 50 1000 600],'Color', 'w');