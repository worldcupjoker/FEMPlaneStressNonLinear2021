% Get real B matrix. 3 x 16
function b = getRealB(B)
b = zeros(3, 16);
for i = 1 : 16
    if (rem(i, 2) == 1)
        b(1, i) = B(1, int16((i + 1) / 2));
        b(2, i) = 0;
        b(3, i) = B(2, int16((i + 1) / 2));
    else
        b(1, i) = 0;
        b(2, i) = B(2, int16(i / 2));
        b(3, i) = B(1, int16(i / 2));
    end
end
return
end