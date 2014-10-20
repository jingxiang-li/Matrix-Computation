%% QR_hh: Decompose A = Q * R using Householder transformation
%% Input:   A, matrix A(m * n) where m > n
%% Output:  A, where upper tri is the R, 
%%          lower tri are the corresponding v in house trans
function [A] = QR_hh(A)
    
    [m, n] = size(A);

    for j = 1 : n
        [v, beta] = house(A(j : m, j));
        tmp = m - j + 1;
        I = eye(tmp);
        A(j:m, j:n) = (I - beta * v * v') * A(j:m, j:n);
        if j < m
            A((j + 1) : m, j) = v(2 : (m - j + 1));
        end
    end

end
