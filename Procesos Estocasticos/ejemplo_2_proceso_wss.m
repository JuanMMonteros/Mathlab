clear all
%close all

% 
Nexp=10000;
Lsim=100;
t=(1:Lsim).';
w0=1;
A=2;

phi = 2*pi*rand(1,Nexp)-pi;
X = A*cos(w0*t+phi);

% Media
media_teo = zeros(size(t));
media_medida = 1/Nexp*sum(X,2);
% figure
% plot(t,media_teo)
% hold all
% plot(t, media_medida,'--')

% Autocorrelacion
% Elijo un t1 y barro t2
t2=25;
t1_v = t;
autocorr_teorica = A^2/2*cos(w0*(t1_v-t2));
aux = X(t2,:).*X;
autocorr_medida = 1/Nexp*sum(aux,2);
figure
plot(t, autocorr_teorica)
hold all
plot(t, autocorr_medida,'--')

% 
