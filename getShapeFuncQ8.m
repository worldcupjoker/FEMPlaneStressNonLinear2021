% Copyright, 2021, (Ted) Zhongkun Zhan, all rights reserved.
% File name: getShapeFuncQ8.m
% Function: Find the shape function (Q8) and it's gradient in
% isoparametric coordination.
% Noding orders: Four corners from the bottom left. Four mid points from
% the bottom.
% Input: [xi] horizontal variable in isoparametric matrix.
%        [eta] vertical variable in isoparametric matrix.
%        [mode] boolean value. true returns the shape function, and false
%        returns the gradient of the shape function.
% Output: 1D or 2D matrix.

function N = getShapeFuncQ8(xi, eta, mode)

% Get the shape function matrix if mode is true.
if mode == true
    N = zeros(1, 8);
    N(1, 1) = -1 / 4 * (1 - xi) * (1 - eta) * (1 + xi + eta);
    N(1, 2) = -1 / 4 * (1 + xi) * (1 - eta) * (1 - xi + eta);
    N(1, 3) = -1 / 4 * (1 + xi) * (1 + eta) * (1 - xi - eta);
    N(1, 4) = -1 / 4 * (1 - xi) * (1 + eta) * (1 + xi - eta);
    N(1, 5) = 1 / 2 * (1 - xi) * (1 + xi) * (1 - eta);
    N(1, 6) = 1 / 2 * (1 + xi) * (1 - eta) * (1 + eta);
    N(1, 7) = 1 / 2 * (1 - xi) * (1 + xi) * (1 + eta);
    N(1, 8) = 1 / 2 * (1 - xi) * (1 - eta) * (1 + eta);
    return
    
% Get the gradient of the shape function, G matrix, if mode is
% false.
else
    N = zeros(2, 8);
    N(1, 1) = -1 / 4 * (-1 + eta) * (2 * xi + eta);
    N(1, 2) = 1 / 4 * (-1 + eta) * (eta - 2 * xi);
    N(1, 3) = 1 / 4 * (1 + eta) * (2 * xi + eta);
    N(1, 4) = -1 / 4 * (1 + eta) * (eta - 2 * xi);
    N(1, 5) = xi * (-1 + eta);
    N(1, 6) = -1 / 2 * (1 + eta) * (-1 + eta);
    N(1, 7) = -xi * (1 + eta);
    N(1, 8) = 1 / 2 * (1 + eta) * (-1 + eta);
    N(2, 1) = -1 / 4 * (-1 + xi) * (xi + 2 * eta);
    N(2, 2) = 1 / 4 * (1 + xi) * (2 * eta - xi);
    N(2, 3) = 1 / 4 * (1 + xi) * (xi + 2 * eta);
    N(2, 4) = -1 / 4 * (-1 + xi) * (2 * eta - xi);
    N(2, 5) = 1 / 2 * (1 + xi) * (-1 + xi);
    N(2, 6) = -eta * (1 + xi);
    N(2, 7) = -1 / 2 * (1 + xi) * (-1 + xi);
    N(2, 8) = eta * (-1 + xi);
    return
end
end