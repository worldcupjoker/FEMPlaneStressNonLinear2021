% Copyright, 2021, (Ted) Zhongkun Zhan, all rights reserved.
% File name: source2PointGaussQ8.m
% N_a: 1 x 8
% N_real_a: 2 x 8

function f_s = source2PointGaussQ8(N_a, N_b, N_c, N_d, detJ_a, detJ_b, detJ_c, detJ_d, L, xy, source, N_real_a, N_real_b, N_real_c, N_real_d)

xa = N_a * xy(:, 1);
ya = N_a * xy(:, 2);
xb = N_b * xy(:, 1);
yb = N_b * xy(:, 2);
xc = N_c * xy(:, 1);
yc = N_c * xy(:, 2);
xd = N_d * xy(:, 1);
yd = N_d * xy(:, 2);

f_s = L.' * (N_real_a.' * detJ_a * source(xa, ya) + N_real_b.' * detJ_b * source(xb, yb) + N_real_c.' * detJ_c * source(xc, yc) + N_real_d.' * detJ_d * source(xd, yd));
return
end