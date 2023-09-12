function [y] = slicer_pam4(x)
%SLICER_PAM4 Summary of this function goes here
%   Detailed explanation goes here
if x<-2
    y=-3;
elseif x<0
    y=-1;
elseif x<2
    y=1;
else
    y=3;
end
end

