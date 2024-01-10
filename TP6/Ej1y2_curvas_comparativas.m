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

% Obtener el título actual y establecer el nuevo tamaño de la letra y formato
title('Respuesta en Frecuencia - Diagrama Bode', 'FontSize', 16, 'FontName', 'Arial');

% Obtener las etiquetas de los ejes y ajustar el tamaño de la letra y formato
xlabel('Frecuency', 'FontSize', 14, 'FontName', 'Helvetica');
ylabel('Phase', 'FontSize', 14, 'FontName', 'Helvetica');

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
bodeopts.PhaseVisible = 'on'; % Oculta la fase
% Obtener la respuesta en frecuencia hasta 7 GHz
bode(H, bodeopts); % Utiliza las opciones personalizadas

% Obtener el título actual y establecer el nuevo tamaño de la letra y formato
title('Respuesta en Frecuencia - Diagrama Bode', 'FontSize', 16, 'FontName', 'Arial');

% Obtener las etiquetas de los ejes y ajustar el tamaño de la letra y formato
xlabel('Frecuency', 'FontSize', 14, 'FontName', 'Helvetica');
ylabel('Phase', 'FontSize', 14, 'FontName', 'Helvetica');

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

hold on;
% Ej2-a
Kp = 0.005;
Ki = 0;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
%figure
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

%figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

hold off;

% Añadir leyenda
legend('Kp = 0.001', 'Kp = 0.005', 'Kp = 0.0005');
title('Comparación de Curvas');
grid on;



%% Ej2-b
Kp = 0.05;
Ki = 0; %Kp/1000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

figure
step(H*z/(z-1))

step(z/(z-1))
% Ajusta el ancho de línea después de generar la gráfica

h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 1); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico}

hold on;

% otro
Kp = 0.005;
Ki = 0;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
%figure
step(H*z/(z-1))

step(z/(z-1))
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','lime'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 3); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico


%otro
Kp = 0.0005;
Ki = 0;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

%figure
step(H*z/(z-1))

step(z/(z-1))
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)

grid on % Muestra la cuadrícula en el gráfico

hold off;

% Añadir leyenda
legend('Kp = 0.05', 'Kp = 0.005', 'Kp = 0.0005');
title('Comparación de Curvas');
grid on;

%% Ej2-c
Kp = 0.05;
Ki = Kp/1000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);


% Configurar opciones del Bode
bodeopts = bodeoptions;
bodeopts.MagVisible = 'on';  % Mostrar solo la parte de magnitud
%bodeopts.PhaseVisible = 'off';  % Ocultar la parte de fase

figure
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

hold on;

Kp = 0.005;
Ki = 0;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

%figure
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

Kp = 0.0005;
Ki = 0;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

%figure
bode(H, bodeopts); % Utiliza las opciones personalizadas
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico



hold off;
% Añadir leyenda

legend('Kp = 0.05', 'Kp = 0.005', 'Kp = 0.0005', 'Location', 'northeast');
title('Comparación de Curvas');
grid on;

%%

% Ej3-1a= variacion de Ki
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

hold on;
Kp = 0.05;
Ki = Kp/1000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
%figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

Kp = 0.05;
Ki = Kp/10000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

%figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

hold off;

% Añadir leyenda
legend('Ki = Kp/100', 'Ki = Kp/1000', 'Ki = Kp/10000');
title('Comparación de Curvas');
grid on;% Añadir leyenda


%%
% Ej3-1a= variacion de kp
Kp = 0.07;
Ki = 0.00007; %Kp/1000

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

hold on;
Kp = 0.01;


H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);
%figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

Kp = 0.005;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

%figure
step(H)
% Ajusta el ancho de línea después de generar la gráfica
h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
grid on % Muestra la cuadrícula en el gráfico

hold off;

% Añadir leyenda
legend('Kp = 0.07', 'Kp = 0.01', 'Kp = 0.005');
title('Comparación de Curvas');
grid on;% Añadir leyenda


%% Ej3-c Variacion de Ki
Kp = 0.05;
Ki = Kp/100;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);


figure;
bode(H); 
hold on;

Kp = 0.05;
Ki = Kp/200;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

bode(H); 

Kp = 0.05;
Ki = Kp/10000;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

bode(H);

hold off;

 h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
 set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
 grid on % Muestra la cuadrícula en el gráfico

 
% Añadir leyenda
legend('Ki = Kp/100', 'Ki = Kp/200', 'Ki = Kp/10000');
title('Comparación de Curvas');
grid on;



%% Ej3-c Variacion de Kp
Kp = 0.07;
Ki = 0,00007;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);


figure;
bode(H); 
hold on;

Kp = 0.01;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

bode(H); 

Kp = 0.005;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

bode(H);

hold off;

 h = findall(gcf,'type','line'); % Obtiene todas las líneas en la figura actual
 set(h, 'LineWidth', 2); % Establece el ancho de línea a 2 (o tu valor deseado)
 grid on % Muestra la cuadrícula en el gráfico

 
% Añadir leyenda
legend('Kp = 0.07', 'Kp = 0.01', 'Kp = 0.005');
title('Comparación de Curvas');
grid on;


