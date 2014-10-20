%% GE-Banded: Gaussian Elimination without pivoting for tridiagonal matrix stored as n * 3
%% Input	A: tridiagonal matrix stored as n * 3
%%			b: n-dimensional vector
%% outputs	x: the solution for system Ax = b

function [x] = GE_Band(A, b)

    n = size(b); %% number of rows
    multiplier = 0; %% multiplier used in each step of GE

    %% for loop for Gaussian Elimination
    for i = 1 : (n - 1)
        multiplier = A(i + 1, 1) / A(i, 2);
        A(i + 1, :) = A(i + 1, :) - [A(i, 2 : 3), 0] * multiplier;
        b(i + 1) = b(i + 1) - b(i) * multiplier;
    endfor
    
    %% Call function to solve the upper triangle linear system
    x = UpTri_Band(A, b);

end