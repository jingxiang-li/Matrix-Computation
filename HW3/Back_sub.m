%%  Back_sub: Backward substitution to solve Ax = b where A is upper-triangular
%%  Input: A, b
%%  Output: x
function [x] = Back_sub(A, b)
    
    n = length(b);

    foo = 0;

    for i = n : -1 : 2
        for j = (i - 1) : -1 : 1
            foo = A(j, i) ./ A(i, i);
            A(j, :) = A(j, :) - foo * A(i, :);
            b(j) = b(j) - foo * b(i);
        end
    end

    x = b ./ diag(A);

end
