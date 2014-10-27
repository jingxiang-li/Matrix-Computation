%% UpTri_Band4: Solve a linear system Ax = b, 
%% where A is an upper triangle matrix stored as n * 4
%% Input:   A, b from QR_Given
%% Output:  x, solution to the linear system

function [x] = UpTri_Band4(A, b)

    n = size(b); %% size of b
    multiplier = 0; %% multiplier used in loop

    for i = (n - 1) : (-1) : 1
        multiplier = A(i, 3) / A(i + 1, 2);
        A(i, :) = A(i, :) - multiplier .* [0, A(i + 1, 1 : 3)];
        b(i) = b(i) - multiplier * b(i + 1);
        if (i > 1)
        	multiplier = A(i - 1, 4) / A(i + 1, 2);
        	A(i - 1, :) = A(i - 1, :) - multiplier .* [0, 0, A(i + 1, 1 : 2)];
        	b(i - 1) = b(i - 1) - multiplier * b(i + 1);
        end	
        
    endfor

    x = b ./ A(:, 2);

end