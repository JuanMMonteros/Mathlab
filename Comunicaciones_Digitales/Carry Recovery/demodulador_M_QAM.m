function [ salida ] = demodulador_M_QAM( entrada, M )
% Toma una señal de entrada con simbolos complejos y genera una secuencia
% de bits que se corresponden con dichos simbolos.

symb = qamdemod(entrada, M);
c = de2bi(symb, 'left-msb');
c = c';
salida = reshape(c,log2(M)*length(symb),1); % bits de salida
