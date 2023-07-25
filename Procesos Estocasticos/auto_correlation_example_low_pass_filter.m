clear all
close all


%%
% Genero un filtro con una determinada respuesta al impulso
hf = rcosine(0.25,1,'normal',0.5,2);
% figure
% plot(hf)
% grid on
% xlabel('Discrete time')
% ylabel('Amplitude')

hf_back = conj(hf(end:-1:1));% Notar que para un filtro real y par, es lo mismo
hf_check = conv(hf_back, hf); 

%%
% Genero ruido blanco
% con media 1 y varianza 1
Nexp = 10e3;
Lsim = 1500; % Vector de muestras bastante mas largo que el filtro
x=1+randn(Lsim,Nexp); % cada columna es un nuevo experimento

% Lo paso por un filtro pasa bajo
y=filter(hf,1,x); % cada columna es el resultado de pasar por el filtro

% save('process_data_white.mat','x');
% save('process_data_colored.mat','y');

%% Calculo autocorrelacion
% Rx(tau) = E{ x(t) x*(t-tau) }
% defino g(tau) = x(t) x*(t-tau)
g_tau = zeros(2*length(y(:,1))-1,Nexp); % Guardo el de g(tau)
my=mean(y(:));
for nexp=1:Nexp
    % 
    %g_tau(:,nexp)= 1/Lsim*conv(y(:,nexp),conj(y(end:-1:1,nexp)));
    g_tau(:,nexp)= xcorr(y(:,nexp)-my,'unbiased'); % equivalente a lo anterior
    
end
corr_result = 1/Nexp * sum(g_tau,2);

%% Calculo media
% mu(t) = E{y(t)}
mu_t = 1/Nexp.* sum(y,2);

figure
plot(mu_t)
grid on

%% Calculo varianza (centralizada)
sigma2_t = 1/Nexp.*sum(abs(y-mean(y(:))).^2,2);

figure
plot(sigma2_t)
grid on

%%
figure
title("Individual contributions of g(\tau)")
plot(g_tau(:,1:10))
grid on
xlabel('Discrete Time')
ylabel('Correlation amplitude')

figure
title("Correlation")
plot(corr_result)
grid on
xlabel('Discrete Time')
ylabel('Correlation amplitude')

%%
[vmax,imax] = max(corr_result);
Ntap = length(hf_check);
corr_cmp = corr_result(imax-(Ntap-1)/2:imax+(Ntap-1)/2);

figure
plot(corr_cmp);
hold all
plot(hf_check);
