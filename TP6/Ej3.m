clear all
close all

fs = 16e9; % GHz
Ts = 1/fs;
z = tf('z', Ts);


L = 100000; 
t = [0:L-1].*Ts;


% Ej3-1a=
Kp = 0.07;
Ki = Kp/1000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

Kp = 0.01;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

Kp = 0.005;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

%%
% Ej3-1b=
Kp = 0.05;
Ki = Kp/100;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

Ki = Kp/1000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

Ki = Kp/10000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

%% Ej3-b
Kp = 0.05;
Ki = Kp/1000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);


step(z/(z-1));
hold on;
step(H*z/(z-1));

% Configuraci�n de la gr�fica
h = findall(gcf,'type','line');
set(h(1),'LineWidth',3,'LineStyle','-','Color','blue'); 

set(h(2),'LineWidth',1,'LineStyle','--','Color','red'); 
legend('Input PLL','Output PLL');
title('Ramp Response');
grid on;

%% Ej3-c
bodeopts = bodeoptions;
bodeopts.FreqUnits = 'Hz'; % Establece las unidades de frecuencia en Hz
bodeopts.MagVisible = 'on'; % Muestra solo la magnitud
bodeopts.PhaseVisible = 'off'; % Oculta la fase
Kp = 0.5;
Ki = Kp/10000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico


Kp = 0.05;
Ki = Kp/1000;

figure
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico


Kp = 0.5;
Ki = Kp/1;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

figure
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico
%%
%EJ3D
Kp = 0.5;
Ki = Kp/10000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H1 = feedback(H0*NCO,z^-0);
figure
pzmap(H, 'x');
hold on;

Kp = 0.5;
Ki = Kp/1;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
pzmap(H, 'x');
