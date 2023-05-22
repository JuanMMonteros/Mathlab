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

P = 4;          % Number of periods for the tone

%% Processing

% Time vector
n_samples = P*OVR;
t_v = (0:n_samples-1)*ts;

% Signals creation
x1_v = sin(2*pi*f0*t_v);
x2_v = cos(2*pi*f0*t_v);
x3_v = x2_v + 1j*x1_v;
x4_v = x2_v - 1j*x1_v;

% Spectrum
NFFT = n_samples;
X1_v = fft(x1_v,NFFT)./n_samples;
X2_v = fft(x2_v,NFFT)./n_samples;
X3_v = fft(x3_v,NFFT)./n_samples;
X4_v = fft(x4_v,NFFT)./n_samples;

% Frequency vector
w_v = (-NFFT/2:NFFT/2-1)*2*pi/NFFT;

%% Plots
fz = 15;

% Signals in time
figure;
plot(t_v/1e-6,x1_v,'-r','Linewidth',2)
hold on;
plot(t_v/1e-6,x2_v,'--b','Linewidth',2)

title('Signals in time','Interpreter','latex','FontSize', fz);
leg_c = {'Sine','Cosine'};
legend(leg_c,'Location','ne','Interpreter','latex','FontSize', fz-2);
xlabel('Time [useg]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
grid on

set(gcf, 'Position', [50 50 400 400],'Color', 'w');

% Signals in freq
figure;

plot(w_v,fftshift(abs(X1_v)),'-r','Linewidth',2)
hold on;
plot(w_v,fftshift(abs(X2_v)),'--b','Linewidth',2)
plot(w_v,fftshift(abs(X3_v)),'--k','Linewidth',1.5)
plot(w_v,fftshift(abs(X4_v)),'-.m','Linewidth',1.5)
title('Signals in frequency','Interpreter','latex','FontSize', fz);
leg_c = {'Sin','Cos','Cos+jSin','Cos-jSin'};
legend(leg_c,'Location','nw','Interpreter','latex','FontSize', fz-2);
xlabel('Angular discrete frequency', 'Interpreter','latex','FontSize', fz);
ylabel('Module', 'Interpreter','latex','FontSize', fz);
grid on; xlim([-pi,pi])

set(gcf, 'Position', [500 50 600 600],'Color', 'w');