clear all
close all

nexp = .1e6;
nsamples = 1e3;

esperanza_por_definicion = zeros(nsamples, 1);
varianza_por_definicion = zeros(nsamples, 1);
figure
for n=1:nexp
    gaussian_noise = randn(nsamples,1);
    wiener_process = cumsum(gaussian_noise);
    if (n<10)
        plot(wiener_process);
        hold all
        
        mean(wiener_process)
    end
    
    esperanza_por_definicion = esperanza_por_definicion+wiener_process/nexp;
    varianza_por_definicion = varianza_por_definicion + 1/nexp*wiener_process.^2;
end

xlabel('Tiempo discreto')
ylabel('Amplitud del ruido')
grid on

plot(esperanza_por_definicion,'--k','LineWidth',2.5)
plot(varianza_por_definicion, '--r','LineWidth',2.5)