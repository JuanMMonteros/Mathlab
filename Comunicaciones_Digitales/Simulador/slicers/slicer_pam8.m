function [y] = slicer_pam8(x)
%SLICER_PAM4 Summary of this function goes here
%   Detailed explanation goes here
if x<-6
    y=-7;
elseif x<-4
    y=-5;
elseif x<-2
    y=-3;
elseif x<0
    y=-1;
elseif x<2
    y=1;
elseif x<4
    y=3;
elseif x<6
    y=5;
else
    y=7;
end
end

