function prod = q_prod(q,p)

p1 = p(1);
p2 = p(2);
p3 = p(3);
p4 = p(4);

% transform quaternion product into matrix product
R = [p4 -p3 p2 p1; p3 p4 -p1 p2; -p2 p1 p4 p3; -p1 -p2 -p3 p4];
prod = R*q;

end