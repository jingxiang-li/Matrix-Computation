%% lanczos Tridiagonalization, Given Matrix, AQ = QT, where Q is orthonormal and T is a tridiagonal Matrix
%% Input: n * n matrix A
%% output: orthonormal Q and tridiagonal T, s.t. AQ = QT

function [Q T] = lanczosTri (A, eps0)
    n = size(A, 1);
    beta = zeros(n + 1 ,1);
    alpha = zeros(n + 1, 1);
    Q = zeros(n, n + 1);
    R = Q;
    v = (rand(n, 1) - 0.5) * 2;
    v = v / norm(v, 2);

    k = 1;
    beta(1) = 1;
    R(:, 1) = v;

    while ((k == 1 || abs(beta(k)) > eps0) && k < n + 1)
        Q(:, k + 1) = R(:, k) / beta(k);
        k = k + 1;
        alpha(k) = Q(:, k)' * A * Q(:, k);
        R(:, k) = (A - alpha(k) * eye(size(A))) * Q(:, k) - beta(k - 1) * Q(:, k - 1);
        beta(k) = norm(R(:, k), 2);
        if abs((R(:, k) / beta(k))' * Q(:, 2)) > eps0
            break;
        end
    end

    Q = Q(:, 2 : (k));
    alpha = alpha(2 : (k));
    beta = beta(2 : (k - 1));
    
    T = zeros(k - 1);
    for i = 1 : k - 1
        T(i, i) = alpha(i);
        if (i < k - 1)
            T(i, i + 1) = beta(i);
            T(i + 1, i) = beta(i);
        end
    end
end