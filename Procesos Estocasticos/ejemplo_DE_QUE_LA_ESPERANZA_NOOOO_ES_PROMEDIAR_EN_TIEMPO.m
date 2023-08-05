clear all
close all

%
Nexp=1000;
fs=100e3;
Ts=1/fs;
f0=2e3;
Tend=5/f0;
t = (0:Ts:Tend-Ts).';

% Genero un ruido iid de varianza 1 y media 0
% randn() me da el resultado de una VA gaussiana de var 1 y media 0
n_t = 1*randn(length(t), Nexp);

x_t = cos(2*pi*f0*t) + n_t;

% Ejemplo, ploteo la ocurrencia del proceso para el experimento 10 y 25
figure
plot(t, x_t(:, 10))
hold all
plot(t, x_t(:, 25))
% Calculo la media = E{X(t)}
media = sum(x_t,2)/Nexp; % Equivalente mean(x_t,2)
plot(t, media, '--k', 'LineWidth',2)
plot(t, cos(2*pi*f0*t), '--r', 'LineWidth',1.2)

% Moraleja: la media puede depender del tiempo
% NO ES UN PROMEDIO A LO LARGO DEL TIEMPO
