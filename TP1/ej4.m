%%
clear all
close all
%%

NPOINTS = 256;
% omega = -pi:2*pi/NPOINTS:pi-2*pi/NPOINTS;
omega = [-1:2/NPOINTS:1-2/NPOINTS].*pi;
n = [0:NPOINTS-1]; % Ojo con el -1!
% definimos alpha
alpha = 0.9;
% Respuesta en frecuencia obtenida analíticamente
H_analitica = (1-alpha)./(1-alpha.*exp(-j.*omega));
% Respuesta al impulso obtenida
h = (1-alpha).*alpha.^n;
% Cálculo de la respuesta en frecuencia utilizando FFT
H_fft = fftshift(fft(h, NPOINTS));

% Gráfico de la respuesta en frecuencia obtenida numéricamente y analíticamente
% Gráfico de la respuesta en frecuencia obtenida numéricamente y analíticamente
figure;
subplot(211); %Amplitud
plot(omega, abs(H_fft), '-b', 'Linewidth', 2);
hold on;
plot(omega, abs(H_analitica), '--r', 'Linewidth', 2);
xlabel('\omega Frequency [rad/s]')
ylabel('Amplitude of H(e^{j\omega})')
legend('FFT', 'Analytical', 'Location', 'best');
grid on;

subplot(212); %Fase
plot(omega, angle(H_fft), '-b', 'Linewidth', 2);
hold on;
plot(omega, angle(H_analitica), '--r', 'Linewidth', 2);
xlabel('\Omega Frequency [rad/s]')
ylabel('Argument of H(e^{j\omega})')
legend('FFT', 'Analytical', 'Location', 'best');
grid on;
%%
alphas = [0.5,0.65,0.8,0.99] % Valores de ? a variar
figure;
hold on;
for i = 1:length(alphas)
    alpha = alphas(i);
    H = (1-alpha)./(1-alpha.*exp(-1j.*omega));
    plot(omega, abs(H), 'LineWidth', 2);
end

xlabel('\omega Frequency [rad/s]');
ylabel('Amplitude of H(e^{j\omega})');
legendCell = cellstr(num2str(alphas', 'Alpha = %.2f'));
legend(legendCell, 'Location', 'best');
grid on;