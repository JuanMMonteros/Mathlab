%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero
%-----------------------------------------------------------------------------%

clear 
close all

%% Basic TX QPSK 
L = 10000;      % Simulation Length
BR = 32e9;      % Baud Rate
N = 4;          % Oversampling rate
rolloff = 0.5;  % Pulse shaping rolloff
h_taps = 101;   % Pulse shaping taps

fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 consecutive samples at Tx output

% Two symbols generation (+1,-1) for QPSK
xi = 2*randi([0,1],L,1)-1; % Real
xq = 2*randi([0,1],L,1)-1; % Imag
x = xi + 1j*xq;

% Upsampling to change sampling rate
xup = upsample(x,N);

% Filter to interpolate the signal
h = raised_cosine(BR/2, fs, rolloff, h_taps, 0);
h_delay = (h_taps-1)/2;
yup = filter(h,1,xup);
%yup = filter(h,1,[xup; zeros(h_taps, 1)]);
%yup = yup(1+h_delay:end);

%% Hasta aca, transmiti la envolvente compleja
% El receptor optimo para este caso es una llave que muestrea a kT segundos
% Esto equivale a decimar xN (pero hay que encontrar la mejor fase de
% muestreo)
for To=0:N-1
    y_rx = yup(To+1:N:50*N);
    scatterplot(y_rx)
    title(sprintf("%d", To));
    set(gcf, 'Position', [50 50 500 500],'Color', 'w');
end
%% Truquito: SI COMPENSO EL DELAY DEL PULSE SHAPING SIEMPRE LA FASE CORRECTA ES LA 0

%%
To = 2;
y_rx = yup(To+1:N:end);
scatterplot(y_rx)
set(gcf, 'Position', [50 50 500 500],'Color', 'w');