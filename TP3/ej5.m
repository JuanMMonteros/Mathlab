Lsim=1000e3;
insd = 1e-9; % nV/sqrt(Hz)
fs = 50e6; %Hz
r = 1; % ohm
No= insd^2/r; % W/Hz
Po= No*fs;
fz=16;
xn = sqrt(Po)*randn(Lsim,1);

NFFT = 1024;
win_v = hanning(128);
noverlap = 0;
[pre,f]=pwelch(xn,hanning(128),0,1024,fs,'twoside');
figure
plot(f,pre,'-k','Linewidth',1.5);
hold on;
grid on
ylim([0,2e-18]);
tit = sprintf('PSDs');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('frequency [Hz]', 'Interpreter','latex','FontSize', fz);
ylabel('PSD amplitude', 'Interpreter','latex','FontSize', fz);
legend({'xn(n)'}, 'Interpreter','latex','FontSize', fz-2);
%% Ejercio B
fs = 100e6;
Po= No*fs;
xn = sqrt(Po)*randn(Lsim,1);
NFFT = 1024;
win_v = hanning(128);
noverlap = 0;
[pre,f]=pwelch(xn,hanning(128),0,1024,fs,'twoside');
figure
plot(f,pre,'-k','Linewidth',1.5);
hold on;
grid on
ylim([0,2e-18]);
tit = sprintf('PSDs');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('frequency [Hz]', 'Interpreter','latex','FontSize', fz);
ylabel('PSD amplitude', 'Interpreter','latex','FontSize', fz);
legend({'xn(n)'}, 'Interpreter','latex','FontSize', fz-2);
