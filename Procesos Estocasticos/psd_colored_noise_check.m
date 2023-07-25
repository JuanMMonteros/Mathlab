clear all
close all

%% Proceso WGN
Lsim=1000e3;
mu_x = 1;
var_x=2;
xn = mu_x + sqrt(var_x)*randn(Lsim,1);

h_lpf = 1/25*ones(25,1);

yn=filter(h_lpf,1,xn);

figure
plot(xn);
hold all
plot(yn,'--');



% Inc. A
NFFT=1024; %'twoside'
[Pxx,f]=pwelch(xn-mu_x,hanning(128),0,1024,'twoside');
figure
plot(f,Pxx)
grid on

% Inc. B
[Pyy,f]=pwelch(yn-mean(yn),hanning(128),0,1024,'twoside');
hold all
plot(f,Pyy)
grid on

potencia_y_psd = sum(Pyy)*(2*pi/NFFT);
potencia_y_var = var(yn);

%%
H_f = fft(h_lpf,NFFT);
S_y_teorico_0 = var_x/2/pi.*abs(H_f).^2;
S_y_teorico_1 = Pxx       .*abs(H_f).^2;

figure
plot(f,Pyy);
hold all
plot(f,S_y_teorico_0);
plot(f,S_y_teorico_1);
grid on
legend("Procesada","Teorica posta", "Teorica hibrida")





