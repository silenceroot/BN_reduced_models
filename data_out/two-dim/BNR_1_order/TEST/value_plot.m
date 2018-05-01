rho  = load('./RHO.csv');
u    = load('./U.csv');
v    = load('./V.csv');
p    = load('./P.csv');
phi_a= load('./PHI_a.csv');
z_a  = load('./Z_a.csv');
phi_b= load('./PHI_b.csv');
z_b  = load('./Z_b.csv');
temper=load('./TEMPER.csv');
d_x=1;
d_y=1;
num_x=size(rho,2);
num_y=size(rho,1);
X = (1:num_x)*d_x;
Y = (1:num_y)*d_y;
[X, Y] = meshgrid(X, Y);
figure(1);
surf(X, Y, rho);
title('rho');
figure(2);
surf(X, Y, p);
title('p');
figure(3);
surf(X, Y, u)
title('u')
figure(4)
surf(X, Y, v);
title('v');
figure(5);
surf(X, Y, phi_a);
title('phi_a');
figure(6);
surf(X, Y, phi_b);
title('phi_b');
figure(7);
surf(X, Y, z_a);
title('z_a');
figure(8);
surf(X, Y, z_b);
title('z_b');
figure(9);
surf(X, Y, temper);
title('temper');
