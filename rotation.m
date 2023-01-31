function A = rotation(q)
% creates the rotation matrix from navigation frame to body frame knowing
% the quaternion q = [q1 q2 q3 q4]; (q4 is the real part)
if length(q) ~= 4
    A = NaN(4);

else
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    
    A = [q1^2-q2^2-q3^3+q4^4 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4);
        2*(q1*q2-q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3+q1*q4);
    2*(q1*q3+q2*q4) 2*(q2*q3-q1*q4) -q1^2-q2^2+q3^2+q4^2]; 
end
end