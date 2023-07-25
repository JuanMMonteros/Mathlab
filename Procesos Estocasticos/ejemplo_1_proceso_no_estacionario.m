clear all
%close all

% 
Nexp=1000;
Lsim=100;
t=(1:Lsim).';

A = randn(1,Nexp);
B = randn(1,Nexp);
X = A + B.*t + t.^2;

% Media
media_teo = t.^2;
media_medida = 1/Nexp*sum(X,2);
% figure
% plot(t,media_teo)
% hold all
% plot(t, media_medida,'--')

% Autocorrelacion
% Elijo un t1 y barro t2
t1=50;
t2_v = t;
autocorr_teorica = 1+t1*t2_v+t1^2*t2_v.^2;
aux = X(t1,:).*X;
autocorr_medida = 1/Nexp*sum(aux,2);
figure
plot(t, autocorr_teorica)
hold all
plot(t, autocorr_medida,'--')


