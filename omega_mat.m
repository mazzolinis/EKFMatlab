function W = omega_mat(w)
% matrix used for quaternion integration
w1 = w(1);
w2 = w(2);
w3 = w(3);

W = [0 w3 -w2 w1;
     -w3 0 w1 w2;
     w2 -w1 0 w3;
     -w1 -w2 -w3 0];

end