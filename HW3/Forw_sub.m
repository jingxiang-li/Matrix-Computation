%%  Forw_sub: Backward substitution to solve Ax = b where A is lower-triangular
%%  Input: A, b
%%  Output: x
function [x] = Forw_sub(A, b)
    
    n = length(b);

    foo = 0;

    for i = 1 : (n - 1)
        for j = (i + 1) : n
            foo = A(j, i) ./ A(i, i);
            A(j, :) = A(j, :) - foo * A(i, :);
            b(j) = b(j) - foo * b(i);
        end
    end

    x = b ./ diag(A);
    
end
