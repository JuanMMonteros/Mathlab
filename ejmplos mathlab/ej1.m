clear all
close all

%% Dada la eq en diff y[n] = alpha*y[n-1] + (1-alpha)*x[n]
% Calcular H(omega) y graficarla
% Paso 1 - Calculo analitico de H(omega)
% H(omega) = (1-alpha)/(1-alpha*e(-j*omega))

%Defino omega
NPOINTS = 256;
% omega = -pi:2*pi/NPOINTS:pi-2*pi/NPOINTS;
omega = [-1:2/NPOINTS:1-2/NPOINTS].*pi;
% .* significa a cada elemento del arreglo, multiplicalo por pi

% Calculo H(omega)
alpha = 0.9;
H = (1-alpha)./(1-alpha.*exp(-j.*omega));
% Con F5 Corro, con F10 avanzo en los bps

% Ploteo H vs omega
figure
subplot 211 % 2 filas, 1 columna, el activo es 1
plot(omega, abs(H), '-b','Linewidth', 2)
xlabel('\omega Frequency [rad/s]')
ylabel('Amplitude of H(e^{j\omega})')
grid on

subplot 212 % 2 filas, 1 columna, el activo es 2
plot(omega, angle(H), '-r','Linewidth', 2)
xlabel('\Omega Frequency [rad/s]')
ylabel('Argument of H(e^{j\omega})')
grid on

%% 
% Paso 2: calcular h(n) en Matlab usando IFFT y comprobar que la rta al
% impulso esta bien
% H(omega) es la FFT de h[n] = (1-alpha)*alpha^k

% Comprobemos
% Como definimos a omega entre [-pi, pi), tenemos que acomodar la rta en
% fcia para que exista entre [0,2pi), usando fftshift()
h0_calc = ifft(fftshift(H));
% Nosotros sabemos que la rta al imp es h[n] = (1-alpha)*alpha^k
n = [0:length(h0_calc)-1]; % Ojo con el -1!
h0_teorica = (1-alpha).*alpha.^n; % El u[n] esta implicito en arrancar el tiempo en 0

figure
plot(n,h0_calc, '--b', 'Linewidth', 2);
hold on
plot(n,h0_teorica, 'ro', 'Markersize', 8);
xlabel('Discrete time')
ylabel('Amplitude of Impulse Response')
grid on
xlim([0,50])
legend('h[n] calc en Matlab', 'h[n] calc analitico')

%% Paso 3: a partir de h[n], calcular en MatLab H(omega) usando FFT
H_matlab = fftshift(fft(h0_teorica));
figure
plot(omega,abs(H_matlab));
xlabel('Frequency \omega')
ylabel('Amplitude of H')
grid on 
% Es igual al modulo de la primera figura