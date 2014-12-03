close all; clear; clc;
cd '~/Documents/Courses/CSCI-5304/CSCI-5304-HW/HW5/code/';

%% 4.a
n = 5
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
iterMax = 1000;
eps0 = 1e-07;
[eigVec eigVal iterNum resArray Rq] = powerMethod(A, iterMax, eps0);
semilogy(1:iterNum, abs(Rq(1:iterNum) - eigVal));

n = 10
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
iterMax = 1000;
eps0 = 1e-07;
[eigVec eigVal iterNum resArray] = powerMethod(A, iterMax, eps0);
semilogy(1:iterNum, abs(Rq(1:iterNum) - eigVal));

n = 40
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
iterMax = 1000;
eps0 = 1e-07;
[eigVec eigVal iterNum resArray Rq] = powerMethod(A, iterMax, eps0);
semilogy(1:iterNum, abs(Rq(1:iterNum) - eigVal));
xlabel('iterNum');ylabel('log(rayleighQuotient - finalValue)');title('log(rayleighQuotient - finalValue) VS iterNum');
fontsize=16;
set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
set([gca; findall(gca, 'Type','line')], 'linewidth', 3);
saveas(1, 'p1.pdf');

%% 4.b
n = 5
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
iterMax = 1000;
eps0 = 1e-07;
[eigVec eigVal iterNum resArray] = QRI(A, iterMax, eps0);
semilogy(1:iterNum, resArray(1:iterNum));

n = 10
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
iterMax = 1000;
eps0 = 1e-07;
[eigVec eigVal iterNum resArray] = QRI(A, iterMax, eps0);
semilogy(1:iterNum, resArray(1:iterNum));

n = 40
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
iterMax = 1000;
eps0 = 1e-07;
[eigVec eigVal iterNum resArray] = QRI(A, iterMax, eps0);
eigVal
%% 4.b modified to force converge to largest eigenvector

n = 5;
iterFix = 3;
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
iterMax = 1000;
eps0 = 1e-07;
[eigVec eigVal iterNum resArray] = QRI_f(A, iterMax, eps0, iterFix);
semilogy(1:iterNum, resArray(1:iterNum));

n = 10;
iterFix = 5;
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
iterMax = 1000;
eps0 = 1e-07;
[eigVec eigVal iterNum resArray] = QRI_f(A, iterMax, eps0, iterFix);
semilogy(1:iterNum, resArray(1:iterNum));

n = 40;
iterFix = 10;
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
iterMax = 1000;
eps0 = 1e-07;
[eigVec eigVal iterNum resArray] = QRI_f(A, iterMax, eps0, iterFix);
semilogy(1:iterNum, resArray(1:iterNum));

%% 4.c
n = 40;
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
iterMax = 1000;
eps0 = 1e-07;
[eigVec eigVal iterNum resArray] = QRAlgo(A, iterMax, eps0);
iterNum
semilogy(1:iterNum, resArray(1:iterNum));

A * eigVec(:, 1) - eigVal(1, 1) * eigVec(:, 1)


%%5
n = 40;
[Q,R] = qr(rand(n,n));
A = Q * diag([1:n-2,2*n-1,2*n]) * Q';
eps0 = 1e-07;
[eigVec eigVal iterNum resArray] = powerMethod(A, iterMax, eps0);
[Q T] = lanczosTri(A, eps0);
res = zeros(size(T, 1), 1);
for i = 1:size(T, 1)
	eigeig = max(eig(T(1:i, 1:i)));
	res(i) = abs(eigeig - eigVal);
end
semilogy(1:size(T,1), res)
