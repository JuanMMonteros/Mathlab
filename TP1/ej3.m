clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window
%% definicion de mi funcion H
%Defino omega
NPOINTS = 256; %resolucion de la w
% omega = -pi:2*pi/NPOINTS:pi-2*pi/NPOINTS;
omega = [-1:2/NPOINTS:1-2/NPOINTS].*pi;
% .* significa a cada elemento del arreglo, multiplicalo por pi

%Definimos nuestra funcion de tranferencia H(w)
H=1/5*(1+exp(-1j.*omega)+exp(-2j.*omega)+exp(-3j.*omega)+exp(-4j.*omega)); 

%plot Modulo de H(w)
figure 
subplot 211 % 2 filas, 1 columna, el activo es 1
plot(omega, abs(H), '-b','Linewidth', 2)
xlabel('\omega Frequency [rad/s]')
ylabel('Amplitude of H(e^{j\omega})')
grid on

%plot Fase de H(w)
subplot 212 % 2 filas, 1 columna, el activo es 2
plot(omega, angle(H), '-r','Linewidth', 2)
xlabel('\Omega Frequency [rad/s]')
ylabel('Argument of H(e^{j\omega})')
grid on

%%
n_taps = 25;
h_v = ones(1, n_taps) / n_taps;% Respuesta al impulso
H1 = fft(h, NPOINTS); % Respuesta en frecuencia utilizando FFT
H1 = fftshift(H1); % Ajustar la posición de los datos utilizando fftshift

figure;
subplot(211);
plot(omega, abs(H1), '-b', 'Linewidth', 2);
hold on;
plot(omega, abs(H), '--r', 'Linewidth', 2);
xlabel('\omega Frequency [rad/s]')
ylabel('Amplitude of H(e^{j\omega})')
legend('vector respusa al impulso', 'original', 'Location', 'best');
grid on;

subplot(212);
plot(omega, angle(H1), '-b', 'Linewidth', 2);
hold on;
plot(omega, angle(H), '--r', 'Linewidth', 2);
xlabel('\Omega Frequency [rad/s]')
ylabel('Argument of H(e^{j\omega})')
legend('vector respusa al impulso', 'original', 'Location', 'best');
grid on;

%%
figure;
hold on;
for M2 = [4, 8, 16, 32]
    h = ones(1, M2) / M2; % Respuesta al impulso con M2 muestras
    H = fft(h, NPOINTS); % Respuesta en frecuencia utilizando FFT
    H = fftshift(H); % Ajustar la posición de los datos utilizando fftshift
    
    plot(omega, abs(H), 'Linewidth', 2);
end
xlim([0,pi])%Analizaremos solamento w>0
xlabel('\omega Frequency [rad/s]')
ylabel('Amplitude of H(e^{j\omega})')
grid on;
legend('M2 = 4', 'M2 = 8', 'M2 = 16', 'M2 = 32', 'Location', 'best');
