function [ y] = slicer_qpsk(x)
%SLICER_QAM 
% M is the modulation format (4,16)

if real(x)>=0
    yr=1;
else
    yr=-1;
end
if imag(x)>=0
    yi = 1;
else
    yi = -1;
end
y = yr+1j*yi;

end

