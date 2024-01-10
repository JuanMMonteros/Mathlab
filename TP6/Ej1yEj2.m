clear all
close all

fs = 16e9; % GHz
Ts = 1/fs;
z = tf('z', Ts);

%%EJ1,
L = 100000; % Simulation Length/////////
t = [0:L-1].*Ts;

% System modelling
Kp = 0.005;
Ki = Kp/100;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);


% Respuesta en frecuencia mostrando solo la magnitud con límites de eje hasta 160 GHz en Hz
figure
bodeopts = bodeoptions;
bodeopts.FreqUnits = 'Hz'; % Establece las unidades de frecuencia en Hz
bodeopts.MagVisible = 'on'; % Muestra solo la magnitud
bodeopts.PhaseVisible = 'on'; % Oculta la fase
% Obtener la respuesta en frecuencia hasta 7 GHz
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

% System modelling
Kp = 0.005;
Ki = Kp/10000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

figure
bodeopts = bodeoptions;
bodeopts.FreqUnits = 'Hz'; % Establece las unidades de frecuencia en Hz
bodeopts.MagVisible = 'on'; % Muestra solo la magnitud
bodeopts.PhaseVisible = 'off'; % Oculta la fase
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

%%
% Ej2-a
Kp = 0.001;
Ki = 0;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico
% Ej2-a
Kp = 0.005;
Ki = 0;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico
Kp = 0.0005;
Ki = 0;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico
%% Ej2-b
Kp = 0.05;
Ki = Kp/1000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
step(H*z/(z-1))
hold all
step(z/(z-1))
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico
% Ej2-a
Kp = 0.005;
Ki = 0;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
step(H*z/(z-1))
hold all
step(z/(z-1))
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

Kp = 0.0005;
Ki = 0;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

figure
step(H*z/(z-1))
hold all
step(z/(z-1))
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

%% Ej2-c
Kp = 0.05;
Ki = Kp/1000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico


Kp = 0.005;
Ki = 0;

%%Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico


Kp = 0.0005;
Ki = 0;

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
figure
step(H*z/(z-1))
hold all
step(z/(z-1))
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

%% Ej3-c    Variacion de Ki
Kp = 0.05;
Ki = Kp/100;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

% Configurar opciones del Bode
bodeopts = bodeoptions;
bodeopts.MagVisible = 'on';  % Mostrar solo la parte de magnitud
bodeopts.PhaseVisible = 'off';  % Ocultar la parte de fase


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

Kp = 0.05;
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



%% Ej3-c    Variacion de Kp
Kp = 0.07;
Ki = 0,00007;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

% Configurar opciones del Bode
bodeopts = bodeoptions;
bodeopts.MagVisible = 'on';  % Mostrar solo la parte de magnitud
bodeopts.PhaseVisible = 'off';  % Ocultar la parte de fase


figure
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

Kp = 0.01;

figure
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

Kp = 0.005;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

figure
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico


