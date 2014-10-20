%% UpTri_Band: Solve a linear system Ax = b, 
%% where A is an upper triangle matrix stored as n * 3 (obtained from GE_Band)
%% Input:   A, b from GE_Band
%% Output:  x, solution to the linear system

function [x] = UpTri_Band(A, b)

    n = size(b); %% size of b
    multiplier = 0; %% multiplier used in loop

    for i = (n - 1) : (-1) : 1
        multiplier = A(i, 3) / A(i + 1, 2);
        A(i, :) = A(i, :) - multiplier .* [0, A(i + 1, 1 : 2)];
        b(i) = b(i) - multiplier * b(i + 1);
    endfor

    x = b ./ A(:, 2);

end