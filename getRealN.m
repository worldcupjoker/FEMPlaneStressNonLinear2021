% Get real N. 2 x 16
function N_real = getRealN(N) % N, 1 x 8.
N_real = zeros(2, 16);
for i = 1 : 16
    if (rem(i, 2) == 1)
        N_real(1, i) = N(1, int16((i + 1) / 2));
        N_real(2, i) = 0;
    else
        N_real(1, i) = 0;
        N_real(2, i) = N(1, int16(i / 2));
    end
end
return
end