rho  = load('./RHO.csv');
u    = load('./U.csv');
v    = load('./V.csv');
p    = load('./P.csv');
phi_a= load('./PHI_a.csv');
z_a  = load('./Z_a.csv');
phi_b= load('./PHI_b.csv');
z_b  = load('./Z_b.csv');
d_x=1;
d_y=1;
num_x=size(rho,2);
num_y=size(rho,1);
X = (1:num_x)*d_x;
Y = (1:num_y)*d_y;
[X, Y] = meshgrid(X, Y);
figure(1);
surf(X, Y, rho(num_x:-1:1,:));
title('rho');
figure(2);
surf(X, Y, p(num_x:-1:1,:));
title('p');
figure(3);
surf(X, Y, u(num_x:-1:1,:))
title('u')
figure(4)
surf(X, Y, v(num_x:-1:1,:));
title('v');
figure(5);
surf(X, Y, phi_a(num_x:-1:1,:));
title('phi_a');
figure(6);
surf(X, Y, phi_b(num_x:-1:1,:));
title('phi_b');
figure(7);
surf(X, Y, z_a(num_x:-1:1,:));
title('z_a');
figure(8);
surf(X, Y, z_b(num_x:-1:1,:));
title('z_b');
