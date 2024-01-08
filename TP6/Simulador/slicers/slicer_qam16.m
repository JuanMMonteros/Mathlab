function [y] = slicer_qam16(x)
%SLICER_QAM16 Summary of this function goes here
%   Detailed explanation goes here
y=slicer_pam4(real(x))+1j.*slicer_pam4(imag(x));
end

