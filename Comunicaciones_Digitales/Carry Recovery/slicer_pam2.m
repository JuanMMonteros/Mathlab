function [y] = slicer_pam2(x)
y=-1*ones(size(x));
mpos = find(x>0);
y(mpos) = 1;

end

