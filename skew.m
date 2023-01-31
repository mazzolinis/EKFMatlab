function  M = skew(v)
% creates the skew-symmetric matrix that represents cross product
if length(v) == 3
    M = [0 -v(3) v(2);
        v(3) 0 -v(1);
        -v(2) v(1) 0];
else
    M = NaN(3);
end
end