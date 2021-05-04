% Copyright, 2021, (Ted) Zhongkun Zhan, all rights reserved.
% File name: force2PointGaussQ8.m
% flux: function handle
% It seems that this function has to be on the path (a seperate file that
% is directly visible to matlab) in order to receive a function handle.
% https://www.mathworks.com/matlabcentral/answers/173155-using-function-handles-for-inputs

function f_n = force2PointGaussQ8(xy, L, recordState, flux)

if recordState == 2
    eta = -1;
    
    % two point gaussian quadrature
    % _xi, _eta are partial derivatives.
    N1 = getShapeFuncQ8(-1 / 3 ^ 0.5, eta, true);
    N1_xi = getShapeFuncQ8(-1 / 3 ^ 0.5, eta, false);
    x1 = N1 * xy(:, 1);
    x1_xi = N1_xi(1, :) * xy(:, 1);
    y1 = N1 * xy(:, 2);
    y1_xi = N1_xi(1, :) * xy(:, 2);
    N2 = getShapeFuncQ8(1 / 3 ^ 0.5, eta, true);
    N2_xi = getShapeFuncQ8(1 / 3 ^ 0.5, eta, false);
    x2 = N2 * xy(:, 1);
    x2_xi = N2_xi(1, :) * xy(:, 1);
    y2 = N2 * xy(:, 2);
    y2_xi = N2_xi(1, :) * xy(:, 2);
    N_real_1 = getRealN(N1);
    N_real_2 = getRealN(N2);
    f_n = L.' * (N_real_1.' * flux(x1, y1) * (x1_xi ^ 2 + y1_xi ^ 2) ^ 0.5 + N_real_2.' * flux(x2, y2) * (x2_xi ^ 2 + y2_xi ^ 2) ^ 0.5);
    return
end

if recordState == 6
    xi = 1;
    
    % two point gaussian quadrature
    N1 = getShapeFuncQ8(xi, -1 / 3 ^ 0.5, true);
    N1_eta = getShapeFuncQ8(xi, -1 / 3 ^ 0.5, false);
    x1 = N1 * xy(:, 1);
    x1_eta = N1_eta(2, :) * xy(:, 1);
    y1 = N1 * xy(:, 2);
    y1_eta = N1_eta(2, :) * xy(:, 2);
    N2 = getShapeFuncQ8(xi, 1 / 3 ^ 0.5, true);
    N2_eta = getShapeFuncQ8(xi, 1 / 3 ^ 0.5, false);
    x2 = N2 * xy(:, 1);
    x2_eta = N2_eta(2, :) * xy(:, 1);
    y2 = N2 * xy(:, 2);
    y2_eta = N2_eta(2, :) * xy(:, 2);
    N_real_1 = getRealN(N1);
    N_real_2 = getRealN(N2);
    f_n = L.' * (N_real_1.' * flux(x1, y1) * (x1_eta ^ 2 + y1_eta ^ 2) ^ 0.5 + N_real_2.' * flux(x2, y2) * (x2_eta ^ 2 + y2_eta ^ 2) ^ 0.5);
    return
end

if recordState == 12
    eta = 1;
    
    % two point gaussian quadrature
    N1 = getShapeFuncQ8(-1 / 3 ^ 0.5, eta, true);
    N1_xi = getShapeFuncQ8(-1 / 3 ^ 0.5, eta, false);
    x1 = N1 * xy(:, 1);
    x1_xi = N1_xi(1, :) * xy(:, 1);
    y1 = N1 * xy(:, 2);
    y1_xi = N1_xi(1, :) * xy(:, 2);
    N2 = getShapeFuncQ8(1 / 3 ^ 0.5, eta, true);
    N2_xi = getShapeFuncQ8(1 / 3 ^ 0.5, eta, false);
    x2 = N2 * xy(:, 1);
    x2_xi = N2_xi(1, :) * xy(:, 1);
    y2 = N2 * xy(:, 2);
    y2_xi = N2_xi(1, :) * xy(:, 2);
    N_real_1 = getRealN(N1);
    N_real_2 = getRealN(N2);
    f_n = L.' * (N_real_1.' * flux(x1, y1) * (x1_xi ^ 2 + y1_xi ^ 2) ^ 0.5 + N_real_2.' * flux(x2, y2) * (x2_xi ^ 2 + y2_xi ^ 2) ^ 0.5);
    return
end

if recordState == 4
    xi = -1;
    
    % two point gaussian quadrature
    N1 = getShapeFuncQ8(xi, -1 / 3 ^ 0.5, true);
    N1_eta = getShapeFuncQ8(xi, -1 / 3 ^ 0.5, false);
    x1 = N1 * xy(:, 1);
    x1_eta = N1_eta(2, :) * xy(:, 1);
    y1 = N1 * xy(:, 2);
    y1_eta = N1_eta(2, :) * xy(:, 2);
    N2 = getShapeFuncQ8(xi, 1 / 3 ^ 0.5, true);
    N2_eta = getShapeFuncQ8(xi, 1 / 3 ^ 0.5, false);
    x2 = N2 * xy(:, 1);
    x2_eta = N2_eta(2, :) * xy(:, 1);
    y2 = N2 * xy(:, 2);
    y2_eta = N2_eta(2, :) * xy(:, 2);
    N_real_1 = getRealN(N1);
    N_real_2 = getRealN(N2);
    f_n = L.' * (N_real_1.' * flux(x1, y1) * (x1_eta ^ 2 + y1_eta ^ 2) ^ 0.5 + N_real_2.' * flux(x2, y2) * (x2_eta ^ 2 + y2_eta ^ 2) ^ 0.5);
    return
end
end