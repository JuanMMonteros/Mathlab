clear all
close all

fs = 16e9; % Frecuencia de muestreo
Ts = 1/fs;
z = tf('z', Ts);

bodeopts = bodeoptions;
bodeopts.FreqUnits = 'Hz'; % Establece las unidades de frecuencia en Hz
bodeopts.MagVisible = 'on'; % Muestra solo la magnitud
bodeopts.PhaseVisible = 'off'; % Oculta la fase


fs = 16e9; % GHz
Ts = 1/fs;
z = tf('z', Ts);

bodeopts = bodeoptions;
bodeopts.FreqUnits = 'Hz'; % Establece las unidades de frecuencia en Hz
bodeopts.MagVisible = 'on'; % Muestra solo la magnitud
bodeopts.PhaseVisible = 'off'; % Oculta la fase

L = [10, 10, 10];
Kp = [0.01, 0.025, 0.05];
Ki = Kp / 1000;

figure
hold on;

for idx = 1:3
    H0 = Kp(idx) + Ki(idx) * z / (z - 1);
    NCO = z / (z - 1);
    H = feedback(H0 * NCO, z^(-L(idx)));
    bode(H, bodeopts); % Utiliza las opciones personalizadas
end

h = findall(gcf, 'type', 'line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

legend('Kp=0.01', 'Kp=0.025', 'K=0.05', 'Location', 'NorthEast');

%% 4b
Kp = 0.05;
Ki = Kp/1000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-10);
figure
pzmap(H, 'x');
hold on;



H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-20);
pzmap(H, 'x');