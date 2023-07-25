clear all
close all

%%
% Ruido real vs complejo
Lsim=1000e3;
Po = 100e-3;
xi = sqrt(Po/2)*randn(Lsim,1);
xq = sqrt(Po/2)*randn(Lsim,1);

x=xi; %Real
x=xi+1j*xq; % Compleja

figure
fs=16e9;
[pxx,f] = pwelch(x,rectwin(1024/2),0,1024,fs);
plot(f,pxx);
grid on
%ylim([0,1])
hold all
ylim([0,7e-12])

%% Flujo mas normal de simulacion
No = 1e-10; %W/Hz 
fs=32e9;
pot_req = No*fs;
x = sqrt(pot_req/2)*randn(100e3,1) + 1j*sqrt(pot_req/2)*randn(100e3,1);

[pxx,f] = pwelch(x,hanning(1024/2),0,1024,fs,'twoside');
plot(f,pxx);
grid on
hold all

