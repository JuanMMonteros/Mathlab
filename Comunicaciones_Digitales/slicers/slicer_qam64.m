function [y] = slicer_qam64(x)
%SLICER_QAM16 Summary of this function goes here
%   Detailed explanation goes here
y=slicer_pam8(real(x))+1j.*slicer_pam8(imag(x));
end

