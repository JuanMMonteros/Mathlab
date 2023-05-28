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

%%

% This script creates an audio signal composed of two tones with noise

% Parameters
fs = 22050;     % Sampling frequency [Hz]

f1=1000;        % Frequency 1 [Hz]
f2=2000;        % Frequency 2 [Hz]

A1 = 1;         % Amplitude 1
A2 = 1;         % Amplitude 2

time_len = 5;   % Time length [s]

% Time vector
t_v = (0:1/fs:time_len);
n_samples = length(t_v);

% Sweep the level of noise
An_v = [10,1,0.1];

for idx = 1:length(An_v)

    An = An_v(idx);         % Noise amplitude

    % Signal creation
    s_v = A1*sin(2*pi*f1*t_v) + A2*sin(2*pi*f2*t_v) + An*randn(1,n_samples); 

    % FFT
    NFFT = 10*n_samples; % oversampling of 10 in the FFT
    S_v = fft(s_v, NFFT)/n_samples; 
    S_db_v = 20*log(abs(S_v)); 
    f_v = (0:NFFT-1)*fs/NFFT;

    % Play signal
    sound(s_v, fs);

    % -- Plots --
    fz = 15;

    % Signals in time
    figure(1); clf;
    plot(t_v,s_v,'-r','Linewidth',2)
    title('Signal in time','Interpreter','latex','FontSize', fz);
    xlabel('Time [s]', 'Interpreter','latex','FontSize', fz);
    ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
    grid on
    set(gcf, 'Position', [50 50 400 400],'Color', 'w');

    % Signals in freq
    figure(2); clf;
    plot((f_v-fs/2)/1e3,fftshift(S_db_v),'-r','Linewidth',2)
    title('Signal in frequency','Interpreter','latex','FontSize', fz);
    xlabel('Frequency [KHz]', 'Interpreter','latex','FontSize', fz);
    ylabel('Module [dB]', 'Interpreter','latex','FontSize', fz);
    grid on; 
    set(gcf, 'Position', [500 50 600 600],'Color', 'w');

    % Wait to audio player
    pause(time_len+1);
    
end