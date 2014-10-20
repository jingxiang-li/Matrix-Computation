%%  LU_nopvt: LU decomposition without partial pivoting using Gaussian Transformation
%%  Input:      A, square matrix recommended
%%  Output:     L U, L is lower-triangular, U is upper-triangular.

function [L U] = LU_nopvt(A)
    
    n = size(A);
    n = n(1);
    
    M = eye(size(A));
    L = eye(size(A));

    for k = 1 : (n - 1)
        gamma = zeros(n, 1);
        for i = (k + 1) : n
            gamma(i) = A(i, k) ./ A(k, k);
        end
        tmp = zeros(n, 1);
        tmp(k) = 1;
        M = eye(n) - gamma * tmp';
        A = M * A;
        L = L * (eye(n) + gamma * tmp');
    end

    U = A;
end