%% 1
x = 1/2 + eps;
y = 1/2;
format hex;
a1 = x^2 - y^2;
a2 = (x + y) * (x - y);
a1
a2

%% 2.1
A = [-0.001, 1.001; 0.001, -0.001];
b = [1; 0];
epsilon = 1e-04;

%% estimate from Theorem
rho = epsilon * cond(A, inf);
error_1 = 2 * epsilon / (1 - rho) * norm(abs(A) * abs(inv(A)), inf);

%% Actual value
A_new = A + epsilon .* abs(A);
b_new = b + epsilon .* abs(b);
y = inv(A_new) * b_new;
x = [1; 1];
norm(x - y, inf) / norm(x, inf)

%% 2.2 -----------------------------------

A = [-0.001, 1.001; 0.001, -0.001];
b = [1; 0];
epsilon = 1e-04;

%% Actual value
A_new = A + epsilon .* abs(A);
b_new = b + epsilon .* abs([0; 1]);
y = inv(A_new) * b_new;
x = [1; 1];
norm(x - y, inf) / norm(x, inf)

%% Estimate

norm(inv(A), inf) * norm(A, inf) / (1 - norm(inv(A), inf) * norm(epsilon * abs(A), inf))... 
* (norm(epsilon * [0;1], inf) / norm(b, inf) + norm(epsilon * abs(A), inf) / norm(A, inf))

%% 3 ------------------------
%% c

A = [1, 2, -1; 2, 1, 1; -1, 1, -2];
B = A + 0.001 * [1; 0; 0] * [0; 0; 1]';
x = [1; -1; -1];
norm(x, 1) / norm(B * x, 1)

%% 4 ---------------------------------------
n = 10;
A = eye(n);
for i = 1 : n
    for j = 1 : n
        if (i - j == 1)
            A(i, j) = -2;
        end
    end
end
E = -(A - eye(n));

result = zeros(n);
for i = 0 : (n - 1)
    result = result + E ^ i;
end

A * result