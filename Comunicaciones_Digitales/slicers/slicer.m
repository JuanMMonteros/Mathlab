function [y] = slicer(x,M)
%SLICER Summary of this function goes here
%   Detailed explanation goes here
L=length(x);
y=zeros(1,L); % Inicializar y como un vector de ceros
for i=1:L % Corregir el bucle para que itere sobre los índices de x
 if M==2     % PAM-2
        y(i) = real(x(i));
        y(i) = sign(y(i));
        
 elseif M==4
    y(i)=slicer_qpsk(x(i));
elseif M==16
    y(i)=slicer_qam16(x(i));
elseif M==64
    y(i)=slicer_qam64(x(i));
else
    assert(0);
end
end

