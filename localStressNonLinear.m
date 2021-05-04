% Copyright, 2021, (Ted) Zhongkun Zhan, all rights reserved.
% File name: localStressNonLinear.m
% Use it for Q8.
% Find stress at 4 Gaussian quadrature points in an element.

function q = localStressNonLinear(xy, N_a, N_b, N_c, N_d, G_a, G_b, G_c, G_d, uLocal, v0, E0, q0)
q = zeros(4, 6);
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

% eij matrix. 2 x 2.
eij_a = zeros(2, 2);
eij_a(1, 1) = eps_a(1, 1) - 1 / 3 * (eps_a(1, 1) + eps_a(2, 1));
eij_a(2, 2) = eps_a(2, 1) - 1 / 3 * (eps_a(1, 1) + eps_a(2, 1));
eij_a(1, 2) = eps_a(3, 1);
eij_a(2, 1) = eij_a(1, 2);

eij_b = zeros(2, 2);
eij_b(1, 1) = eps_b(1, 1) - 1 / 3 * (eps_b(1, 1) + eps_b(2, 1));
eij_b(2, 2) = eps_b(2, 1) - 1 / 3 * (eps_b(1, 1) + eps_b(2, 1));
eij_b(1, 2) = eps_b(3, 1);
eij_b(2, 1) = eij_b(1, 2);

eij_c = zeros(2, 2);
eij_c(1, 1) = eps_c(1, 1) - 1 / 3 * (eps_c(1, 1) + eps_c(2, 1));
eij_c(2, 2) = eps_c(2, 1) - 1 / 3 * (eps_c(1, 1) + eps_c(2, 1));
eij_c(1, 2) = eps_c(3, 1);
eij_c(2, 1) = eij_c(1, 2);

eij_d = zeros(2, 2);
eij_d(1, 1) = eps_d(1, 1) - 1 / 3 * (eps_d(1, 1) + eps_d(2, 1));
eij_d(2, 2) = eps_d(2, 1) - 1 / 3 * (eps_d(1, 1) + eps_d(2, 1));
eij_d(1, 2) = eps_d(3, 1);
eij_d(2, 1) = eij_d(1, 2);

% Scalar.
epse_a = (2 / 3 * ((eij_a(1, 1)) ^ 2 + (eij_a(2, 2)) ^ 2 + (eij_a(1, 2)) ^ 2 + (eij_a(2, 1)) ^ 2)) ^ 0.5;
epse_b = (2 / 3 * ((eij_b(1, 1)) ^ 2 + (eij_b(2, 2)) ^ 2 + (eij_b(1, 2)) ^ 2 + (eij_b(2, 1)) ^ 2)) ^ 0.5;
epse_c = (2 / 3 * ((eij_c(1, 1)) ^ 2 + (eij_c(2, 2)) ^ 2 + (eij_c(1, 2)) ^ 2 + (eij_c(2, 1)) ^ 2)) ^ 0.5;
epse_d = (2 / 3 * ((eij_d(1, 1)) ^ 2 + (eij_d(2, 2)) ^ 2 + (eij_d(1, 2)) ^ 2 + (eij_d(2, 1)) ^ 2)) ^ 0.5;

v = v0;
E_a = E0 * (1 + epse_a ^ q0);
E_b = E0 * (1 + epse_b ^ q0);
E_c = E0 * (1 + epse_c ^ q0);
E_d = E0 * (1 + epse_d ^ q0);

%%%
% Dr matrix
Dr_a = getDMatrix(v, E_a);
Dr_b = getDMatrix(v, E_b);
Dr_c = getDMatrix(v, E_c);
Dr_d = getDMatrix(v, E_d);

sigma_a = Dr_a * eps_a;
sigma_b = Dr_b * eps_b;
sigma_c = Dr_c * eps_c;
sigma_d = Dr_d * eps_d;

sigma_e_a = 3 * E_a / 2 / (1 + v) * epse_a;
sigma_e_b = 3 * E_b / 2 / (1 + v) * epse_b;
sigma_e_c = 3 * E_c / 2 / (1 + v) * epse_c;
sigma_e_d = 3 * E_d / 2 / (1 + v) * epse_d;

% For testing purpose.
% E_a / (1 - 2 * v0) * (eps_a(1, 1) + eps_a(2, 1));
% E_b / (1 - 2 * v0) * (eps_b(1, 1) + eps_b(2, 1));
% E_c / (1 - 2 * v0) * (eps_c(1, 1) + eps_c(2, 1));
% E_d / (1 - 2 * v0) * (eps_d(1, 1) + eps_d(2, 1));
%%%

% 3 x 1 matrix.
% sigma_a = E_a / (1 + v) * eps_a + v * E_a / (1 + v) / (1 - 2 * v) * (eps_a(1, 1) + eps_a(2, 1)) * [1; 1; 0];
% sigma_b = E_b / (1 + v) * eps_b + v * E_b / (1 + v) / (1 - 2 * v) * (eps_b(1, 1) + eps_b(2, 1)) * [1; 1; 0];
% sigma_c = E_c / (1 + v) * eps_c + v * E_c / (1 + v) / (1 - 2 * v) * (eps_c(1, 1) + eps_c(2, 1)) * [1; 1; 0];
% sigma_d = E_d / (1 + v) * eps_d + v * E_d / (1 + v) / (1 - 2 * v) * (eps_d(1, 1) + eps_d(2, 1)) * [1; 1; 0];

% Store local stress array.
q(1, 3) = sigma_a(1, 1);
q(1, 4) = sigma_a(2, 1);
q(1, 5) = sigma_a(3, 1);
q(1, 6) = sigma_e_a;
q(2, 3) = sigma_b(1, 1);
q(2, 4) = sigma_b(2, 1);
q(2, 5) = sigma_b(3, 1);
q(2, 6) = sigma_e_b;
q(3, 3) = sigma_c(1, 1);
q(3, 4) = sigma_c(2, 1);
q(3, 5) = sigma_c(3, 1);
q(3, 6) = sigma_e_c;
q(4, 3) = sigma_d(1, 1);
q(4, 4) = sigma_d(2, 1);
q(4, 5) = sigma_d(3, 1);
q(4, 6) = sigma_e_d;
return
end