clear all
close all

rand('seed',1); % Setear la semilla es una buena practica

%% Modelo del dado
Nexp = 1000e3;
% rand genera numeros aleatorios UNIFORMES para variable continua
% randi genera numeros aleatorios UNIFORMES para variable discreta

outcomes=zeros(Nexp,1);
for n=1:Nexp
    outcome = randi([1,6],1,1);
    outcomes(n)=outcome;
end
% Todo esto es equivalente a randi([1 6], Nexp,1)

% PMF
[histograma, bins] = hist(outcomes,6);
bins=1:6;
figure;
stem(bins, histograma)
xlim([0,7])
ylim([0,Nexp/6+1000])
grid on
xlabel('Dice outcome')
ylabel('Ocurrences')


pmf = histograma/Nexp;
figure;
stem(bins, pmf)
xlim([0,7])
ylim([0,0.3])
grid on
xlabel('Dice outcome')
ylabel('PMF')


%% Modelo de uniforme continua entre -2 y 5
Nexp = 1000e3;
% rand genera numeros aleatorios UNIFORMES para variable continua
% randi genera numeros aleatorios UNIFORMES para variable discreta

outcomes=zeros(Nexp,1);
for n=1:Nexp
    outcome = 7*rand()-2;
    outcomes(n)=outcome;
end
% Esto es equivalente a rand(Nexp,1)*7-2

% PMF
[histograma, bins] = hist(outcomes,1000);
figure;
plot(bins, histograma,'-o')
xlim([-10,10])
%ylim([0,Nexp/6+1000])
grid on
xlabel('VA outcome')
ylabel('Ocurrences')


step=bins(2)-bins(1);
pdf = histograma/Nexp/step;
figure;
plot(bins, pdf,'-o')
xlim([-10,10])
%ylim([0,0.3])
grid on
xlabel('Dice outcome')
ylabel('PMF')


%% Modelo de Normal media 2, varianza 1.5
Nexp = 1000e3;
% randn --> Normal media 0, varianza=1
gmean=2;
gvar=1.5;
outcomes=zeros(Nexp,1);
for n=1:Nexp
    outcome = sqrt(gvar)*randn()+gmean; % sqrt(VAR) * X + mu, donde N~N(0,1)
    outcomes(n)=outcome;
end
% Es equivalente a sqrt(gvar)*randn(Nexp,1)+gmean

% PMF
[histograma, bins] = hist(outcomes,50);
figure;
plot(bins, histograma,'-o')
xlim([-10,10])
%ylim([0,Nexp/6+1000])
grid on
xlabel('Dice outcome')
ylabel('Ocurrences')


step=bins(2)-bins(1);
pdf = histograma/Nexp/step;
figure;
plot(bins, pdf,'-o')
xlim([-10,10])
%ylim([0,0.3])
grid on
xlabel('VA outcome')
ylabel('PMF')

hold all
xg=bins;
plot(xg, 1/sqrt(2*pi*gvar).*exp(-(xg-gmean).^2/2/gvar),'LineWidth',2.5)
legend("Experiment","Theoretical fit")
    