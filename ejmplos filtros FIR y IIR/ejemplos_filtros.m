clear all
close all


%% requiero convolucionar con una respuesta al impulso h[n] = [0.1,0.2,0.3,0.4,0.3,0.2,0.1];
h = [0.1,0.2,0.3,0.4,0.3,0.2,0.1];
h=h/sum(h); % Forzar que tenga ganancia de DC unitaria
x = randn(100e3,1)+0.8;
% Es tipico es trabajar column-wise
Ldata = length(x);

% Ejemplo filter
tic
y_fir = filter(h, 1, x); % Filtro FIR
toc

% IIR
alpha = 1e-1;
b = [alpha];
a = [1, -(1-alpha)];
y_iir = filter(b,a,x);

figure
plot(x, '--b')
grid on
hold all
plot(y_fir, '-r', 'LineWidth',2.5)
plot(y_iir, '-c', 'LineWidth',1.5)
xlabel('Tiempo Discreto')
ylabel('Amplitud')
legend(["Entrada", "Salida FIR", "Salida IIR"])
plot(0.8*ones(length(y_iir),1), '--k')

%% Averiguar la respuesta en frecuencia de un sistema LTI
[H_iir, f_iir] = freqz(b,a, 64*1024, 1);
[H_fir, f_fir] = freqz(h,1, 64*1024, 1);
%[H, f] = freqz(h,1, 64*1024, 1);
figure
subplot 211
plot(2*pi*f_iir, (abs(H_iir)))
hold all
plot(2*pi*f_fir, (abs(H_fir)))
grid on
xlabel("Discrete Angular Frequency")
ylabel("Module")

subplot 212
plot(2*pi*f_iir, (angle(H_iir)))
hold all
plot(2*pi*f_fir, (angle(H_fir)))
grid on
xlabel("Discrete Angular Frequency")
ylabel("Phase")
% RESUMEN 
% conv
% filter 
% freqz 

%% Ejemplo del filtro FIR
% Filtro de promedios moviles (MAF)
M = 128*16;
h = 1/M.*ones(M,1);
% figure
% stem(h)

y = filter(h,1,x);

figure
plot(x, '--b')
grid on
hold all
plot(y, '-r', 'LineWidth',2.5)
xlabel('Tiempo Discreto')
ylabel('Amplitud')
legend(["Entrada", "Salida"])
plot(0.8*ones(length(y),1), '--k')

%%
M_list = [2,4,8,16,32]*16;

figure
labels={};
for idx=1:length(M_list)
   M = M_list(idx);
   h = 1/M.*ones(M,1);
   
   [H_fir,f_fir] = freqz(h, 1, 8*1024,1);
   plot(2*pi*f_fir, 20*log10(abs(H_fir)), 'LineWidth', 2);
   labels{end+1} = sprintf("M=%d", M);
   hold all
      
end

grid on
xlabel("Frecuencia discreta")
ylabel("Magnitud [dB]")
xlim([0,pi])

plot(2*pi*f_fir, -3*ones(size(f_fir)), '--k')
labels{end+1} = sprintf("-3dB", M);

legend(labels)

%% Evaluacion de una respuesta en frecuencia
% Ejemplo 2.155
alpha = 0.5;
NFFT = 8*1024;
delta_omega = 2*pi/NFFT;
omega = 0:delta_omega:2*pi-delta_omega;
%omega = -pi:delta_omega:pi-delta_omega; % MAAAAl
X = 1 ./ (1- alpha * exp(-1j .*omega));

% Asi lo entiende Matlab
figure
plot(omega, (abs(X)))
xlabel('Freq. Discreta')
ylabel('Margnitud')
grid on

% Asi lo podemos ver nosotros
figure
plot(omega-pi, fftshift(abs(X)))
xlabel('Freq. Discreta')
ylabel('Margnitud')
grid on





