function [y] = slicer(x,M)
%SLICER Summary of this function goes here
%   Detailed explanation goes here
 if M==2     % PAM-2
        y = real(x);
        y(real(y)<0) = -1;
       y(real(y)>0) = 1;
        
 elseif M==4
    y=slicer_qpsk(x);
elseif M==16
    y=slicer_qam16(x);
elseif M==64
    y=slicer_qam64(x);
else
    assert(0);
end
end

