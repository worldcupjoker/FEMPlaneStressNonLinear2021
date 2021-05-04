% Get D matrix.
function D = getDMatrix(v, E)
D = zeros(3, 3);
D(1, 1) = 1;
D(2, 2) = 1;
D(3, 3) = (1 - v) / 2;
D(1, 2) = v;
D(2, 1) = v;
D = D * (E / (1 - v ^ 2));
return
end