% Parámetros
alpha = 0.95;
M = 40;

% Filtro IIR
NPOINTS = 256;
omega = [-1:2/NPOINTS:1-2/NPOINTS].*pi;
H_IIR = (1-alpha)./(1-alpha.*exp(-1j.*omega));

% Filtro FIR
h = ones(1, M)/M;
H_FIR = fftshift(fft(h, NPOINTS));

% Supongamos que tienes la respuesta en frecuencia del filtro IIR almacenada en el vector H_IIR

% Nivel de referencia para calcular el ancho de banda (por ejemplo, -3 dB)
nivelReferencia = 0.707;

% Encontrar las frecuencias de corte donde la respuesta en frecuencia cae al nivel de referencia
indicesCorte = find(abs(H_IIR) <= nivelReferencia);
frecuenciaCorte1 = omega(indicesCorte(1));
frecuenciaCorte2 = omega(indicesCorte(end));

% Calcular el ancho de banda
anchoBanda = abs(frecuenciaCorte2 - frecuenciaCorte1);
fprintf('El ancho de banda aproximado del filtro IIR es: %.2f rad/s\n', anchoBanda);
% Gráficos
figure;
subplot(211);
plot(omega, abs(H_IIR), '-b', 'Linewidth', 2);
hold on;
plot(omega, abs(H_FIR), '--r', 'Linewidth', 2);
xlabel('\omega Frequency [rad/s]');
ylabel('Amplitude');
legend('IIR Filter', 'FIR Filter');
grid on;

