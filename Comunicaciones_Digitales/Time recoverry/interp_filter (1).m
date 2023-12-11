function [ipr] = interp_filter(Ntaps, t0)
rolloff = 0.05 +0.01; % Para evitar los puntos donde el RC no existe
fs=1;fc=1;
Ts=1/fs;
T = 1/fc;

if mod(Ntaps,2)==0
    Ntaps=Ntaps+1; %Fuerzo cant de taps impar
end

t= [-(Ntaps-1)/2:1:(Ntaps-1)/2].*Ts+t0;
t_norm = t./T;
ipr = sinc(t_norm).*( cos(pi.*rolloff.*t_norm) ) ./ (1- (2*rolloff.*t_norm).^2 ); 

end

