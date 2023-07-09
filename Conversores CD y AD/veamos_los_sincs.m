clear all
close all

delays = [-0.8,-0.5,0.25,0,[-0.8,-0.5,0.25]*-1];
delays = [0,0.5];

labels={};
figure
for delay=delays
    %delay=0.5; % esto representa una fraccion de muestra
    % Genero el filtro de la eq 4.65 del libro
    ntaps=13; % Mantenerlo impar para simplificar
    group_delay = (ntaps-1)/2; % Este delay hay que compensarlo porque representa la cola de la convolucion
    nline = -(ntaps-1)/2:(ntaps-1)/2;
    h1 = sinc(nline-delay);
    %h1= my_rcosine(1,1,0.5,ntaps,delay);

    stem(nline, h1);
    hold all
    labels{end+1} = sprintf("D: %2.2f",delay);
end
legend(labels)
    
grid on
xlabel('Discrete time')
title('Interpolation filter')
    