function [s, J] = lkin(k,alpha)
% calculates the direct kinematic and Jacobian of every leg (calculated
% from Wolfram Mathematica
q1 = alpha(1);
q2 = alpha(2);
q3 = alpha(3);

% parameters
x0 = 0.25; % half length of robot back
y0 = 0.15; % half width
a2 = 0.4;
a3 = 0.5;


if k == 1
    s = [x0+(-1).*(a2.*cos(q2)+a3.*cos(q2+q3)).*sin(q1),y0+cos(q1).*(a2.* ...
        cos(q2)+a3.*cos(q2+q3)),a2.*sin(q2)+a3.*sin(q2+q3)];

    J = [(-1).*cos(q1).*(a2.*cos(q2)+a3.*cos(q2+q3)),(-1).*sin(q1).*((-1) ...
        .*a2.*sin(q2)+(-1).*a3.*sin(q2+q3)),a3.*sin(q1).*sin(q2+q3);(-1).* ...
        (a2.*cos(q2)+a3.*cos(q2+q3)).*sin(q1),cos(q1).*((-1).*a2.*sin(q2)+ ...
        (-1).*a3.*sin(q2+q3)),(-1).*a3.*cos(q1).*sin(q2+q3);0,a2.*cos(q2)+ ...
        a3.*cos(q2+q3),a3.*cos(q2+q3)];

elseif k == 2
    s = [x0+(a2.*cos(q2)+a3.*cos(q2+q3)).*sin(q1),(-1).*y0+(-1).*cos(q1).* ...
        (a2.*cos(q2)+a3.*cos(q2+q3)),a2.*sin(q2)+a3.*sin(q2+q3)];

    J = [cos(q1).*(a2.*cos(q2)+a3.*cos(q2+q3)),sin(q1).*((-1).*a2.*sin(q2) ...
        +(-1).*a3.*sin(q2+q3)),(-1).*a3.*sin(q1).*sin(q2+q3);(a2.*cos(q2)+ ...
        a3.*cos(q2+q3)).*sin(q1),(-1).*cos(q1).*((-1).*a2.*sin(q2)+(-1).* ...
        a3.*sin(q2+q3)),a3.*cos(q1).*sin(q2+q3);0,a2.*cos(q2)+a3.*cos(q2+ ...
        q3),a3.*cos(q2+q3)];

elseif k == 3
    s = [(-1).*x0+(a2.*cos(q2)+a3.*cos(q2+q3)).*sin(q1),(-1).*y0+(-1).* ...
        cos(q1).*(a2.*cos(q2)+a3.*cos(q2+q3)),a2.*sin(q2)+a3.*sin(q2+q3)];

    J = [cos(q1).*(a2.*cos(q2)+a3.*cos(q2+q3)),sin(q1).*((-1).*a2.*sin(q2) ...
        +(-1).*a3.*sin(q2+q3)),(-1).*a3.*sin(q1).*sin(q2+q3);(a2.*cos(q2)+ ...
        a3.*cos(q2+q3)).*sin(q1),(-1).*cos(q1).*((-1).*a2.*sin(q2)+(-1).* ...
        a3.*sin(q2+q3)),a3.*cos(q1).*sin(q2+q3);0,a2.*cos(q2)+a3.*cos(q2+ ...
        q3),a3.*cos(q2+q3)];

elseif k == 4
    s = [(-1).*x0+(-1).*(a2.*cos(q2)+a3.*cos(q2+q3)).*sin(q1),y0+cos(q1).* ...
        (a2.*cos(q2)+a3.*cos(q2+q3)),a2.*sin(q2)+a3.*sin(q2+q3)];

    J = [(-1).*cos(q1).*(a2.*cos(q2)+a3.*cos(q2+q3)),(-1).*sin(q1).*((-1) ...
        .*a2.*sin(q2)+(-1).*a3.*sin(q2+q3)),a3.*sin(q1).*sin(q2+q3);(-1).* ...
        (a2.*cos(q2)+a3.*cos(q2+q3)).*sin(q1),cos(q1).*((-1).*a2.*sin(q2)+ ...
        (-1).*a3.*sin(q2+q3)),(-1).*a3.*cos(q1).*sin(q2+q3);0,a2.*cos(q2)+ ...
        a3.*cos(q2+q3),a3.*cos(q2+q3)];

else 
    s = NaN(1,3);
    J = NaN(3);
end
s = s';
end
