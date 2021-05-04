% Copyright, 2021, (Ted) Zhongkun Zhan, all rights reserved.
% File name: localStressLinear.m
% Use it for Q8.
% Find stress at 4 Gaussian quadrature points in an element.

function q = localStressLinear(xy, N_a, N_b, N_c, N_d, G_a, G_b, G_c, G_d, uLocal, v0, E0, q0)
q = zeros(4, 5);
u_n = zeros(16, 1);

for i = 1 : 8
    u_n(i * 2 - 1, 1) = uLocal(i, 1);
    u_n(i * 2, 1) = uLocal(i, 2);
end

% Get x, y coordinates for each quadrature point.
xa = N_a * xy(:, 1);
ya = N_a * xy(:, 2);
q(1, 1) = xa;
q(1, 2) = ya;
xb = N_b * xy(:, 1);
yb = N_b * xy(:, 2);
q(2, 1) = xb;
q(2, 2) = yb;
xc = N_c * xy(:, 1);
yc = N_c * xy(:, 2);
q(3, 1) = xc;
q(3, 2) = yc;
xd = N_d * xy(:, 1);
yd = N_d * xy(:, 2);
q(4, 1) = xd;
q(4, 2) = yd;

% Find the jacobian matrix at each quadrature point.
% Find the flux at each quadrature point.
J_a = G_a * xy;
J_b = G_b * xy;
J_c = G_c * xy;
J_d = G_d * xy;

% Find B matrices. 2 x 8
B_a = inv(J_a) * G_a;
B_b = inv(J_b) * G_b;
B_c = inv(J_c) * G_c;
B_d = inv(J_d) * G_d;

% Real B matrices. 3 x 16
B_real_a = getRealB(B_a);
B_real_b = getRealB(B_b);
B_real_c = getRealB(B_c);
B_real_d = getRealB(B_d);

% Epsilon matrix: eps_xx, eps_yy, gamma_xy
eps_a = B_real_a * u_n;
eps_b = B_real_b * u_n;
eps_c = B_real_c * u_n;
eps_d = B_real_d * u_n;

% D matrix
D = getDMatrix(v0, E0);

% 3 x 1 matrix.
sigma_a = D * eps_a;
sigma_b = D * eps_b;
sigma_c = D * eps_c;
sigma_d = D * eps_d;

% Store local stress array.
q(1, 3) = sigma_a(1, 1);
q(1, 4) = sigma_a(2, 1);
q(1, 5) = sigma_a(3, 1);
q(2, 3) = sigma_b(1, 1);
q(2, 4) = sigma_b(2, 1);
q(2, 5) = sigma_b(3, 1);
q(3, 3) = sigma_c(1, 1);
q(3, 4) = sigma_c(2, 1);
q(3, 5) = sigma_c(3, 1);
q(4, 3) = sigma_d(1, 1);
q(4, 4) = sigma_d(2, 1);
q(4, 5) = sigma_d(3, 1);
return
end