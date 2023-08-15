clc; close all; clear
%%
Lsim=1e9;
Po = 100e-3;
OVR=8;
fs=8e9;
fz=15;

xn = sqrt(Po)*randn(Lsim,1);
%% PDS incizo A 
NFFT = 1024;
win_v = hanning(128);
noverlap = 0;
[pre,f]=pwelch(xn,hanning(128),0,1024,fs,'twoside');
figure
plot(f,pre,'-k','Linewidth',1.5);
hold on;
grid on
tit = sprintf('PSDs');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('frequency [Hz]', 'Interpreter','latex','FontSize', fz);
ylabel('PSD amplitude', 'Interpreter','latex','FontSize', fz);
legend({'xn(n con fs= 8Ghz)'}, 'Interpreter','latex','FontSize', fz-2);
ylim([0,2e-8])
% la densisdad de potenica de ruido es igual a la 
%% PDS incizo  B 
fs=16e6;
[pxx,f]=pwelch(xn,hanning(128),0,1024,fs,'twoside');
figure
plot(f,pxx,'-k','Linewidth',1.5);
hold on;
grid on
tit = sprintf('PSDs');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('frequency [Hz]', 'Interpreter','latex','FontSize', fz);
ylabel('PSD amplitude', 'Interpreter','latex','FontSize', fz);
legend({'xn(n con fs= 16Ghz)'}, 'Interpreter','latex','FontSize', fz-2);
ylim([0,1e-8])
%% incizo C
xi = sqrt(Po/2)*randn(Lsim,1);
xq = sqrt(Po/2)*randn(Lsim,1);
xncmp =xi+1j*xq; % Compleja
[pcmp,f]=pwelch(xncmp,hanning(128),0,1024,fs,'twoside');
figure
plot(f,pcmp,'-k','Linewidth',1.5);
hold on;
plot(f,pxx,'--b','Linewidth',2);
grid on
ylim([0,1e-8])
tit = sprintf('PSDs');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('frequency [Hz]', 'Interpreter','latex','FontSize', fz);
ylabel('PSD amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Ruido Complejo','Ruidos Real'}, 'Interpreter','latex','FontSize', fz-2);


