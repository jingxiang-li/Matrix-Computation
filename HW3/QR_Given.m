%% QR_Given: Solve Ax = b, where A is a tridiagonal matrix stored as n * 3
%% Using Given Rotation
%% Input:   A, n * 3
%%          b, n * 1
%% Output:  A and b, s.t. A * x = b, where A is uppertri stored as n * 3
%%          A is R and b is Q'*b
%%          Q is stored as [c, s, i, j]

function [A b Q] = QR_Given(A, b)

    n = size(b)(1);
    A = [A, zeros(n, 1)];
    A((n - 1) : n, 4) = NaN;
    Q = zeros(n - 1, 4);

    for i = 1 : (n - 1)
        j = i + 1;
        xi = A(i, 2);
        xj = A(j, 1);
        c = xi / sqrt(xi^2 + xj^2);
        s = -xj / sqrt(xi^2 + xj^2);
        A(i, 2) = c * xi - s * xj;
        A(j, 1) = s * xi + c * xj;
        [tmp1 tmp2] = Given_trans(A(i, 3), A(j, 2), c, s);
        A(i, 3) = tmp1;
        A(j, 2) = tmp2;
        if (j < n)
            [tmp1 tmp2] = Given_trans(A(i, 4), A(j, 3), c, s);
            A(i, 4) = tmp1;
            A(j, 3) = tmp2;
        end
        bi = b(i);
        bj = b(j);
        b(i) = c * bi - s * bj;
        b(j) = s * bi + c * bj;
        Q(i, :) = [c, s, i, j];
    end

end