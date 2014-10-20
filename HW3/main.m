close all; clear; clc;
cd /home/jingxiang-li/Documents/Courses/CSCI-5304/CSCI-5304-HW/HW3

m = 10;
A_trid = [ [nan; ones(2*m,1)*m ], (m+1+(-m:m)') , [ones(2*m,1)*m;nan]];
A_full = diag(A_trid(:,2))+diag(A_trid(2:end,1),-1)+diag(A_trid(1:end-1,3),1);
b = ones(size(A_trid, 1), 1);

[GE_Band(A_trid, b) inv(A_full) * b]

%%%%%%%%%%%%%%%%%%%

m_array = 500 : 500 : 5000;
result_time = zeros(size(m_array)(2), 1);
result_res = zeros(size(m_array)(2), 1);
for i = 1 : size(m_array)(2)
    m = m_array(i);
    A_trid = [ [nan; ones(2*m,1)*m ], (m+1+(-m:m)') , [ones(2*m,1)*m;nan]];
    A_full = diag(A_trid(:,2))+diag(A_trid(2:end,1),-1)+diag(A_trid(1:end-1,3),1);
    b = ones(size(A_trid, 1), 1);
    tic();
        x1 = GE_Band(A_trid, b);
        % x1 = A_full \ b;
    result_time(i) = toc();
    result_res(i) = norm(b - A_full * x1, 1);
endfor

[result_time, result_res]

plot(m_array, result_time)
% plot(m_array, result_res)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc

m = 5;
A_trid = [ [nan; ones(2*m,1)*m ], (m+1+(-m:m)') , [ones(2*m,1)*m;nan]];
A_full = diag(A_trid(:,2))+diag(A_trid(2:end,1),-1)+diag(A_trid(1:end-1,3),1);
b = ones(size(A_trid, 1), 1);

[R b1] = QR_Given(A_trid, b);
x = UpTri_Band4(R, b1);
[A_full \ b x]



m = 5;
A_trid = [ [nan; ones(2*m,1)*m ], (m+1+(-m:m)') , [ones(2*m,1)*m;nan]];
A_full = diag(A_trid(:,2))+diag(A_trid(2:end,1),-1)+diag(A_trid(1:end-1,3),1);
b = ones(size(A_trid, 1), 1);

[R b1 Q] = QR_Given_Q(A_trid, b);













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [1, -1; 1, -1.00001];
B = [1, -1; -1, 1.00001];
eig_A = eig(A);
eig_B = eig(B);
max(abs(eig_A)) / min(abs(eig_A))
% > 1.0032
max(abs(eig_B)) / min(abs(eig_B))
% > 4.0000e+05

max(svd(A)) / min(svd(A))
% > 4.0000e+05
max(svd(B)) / min(svd(B))
% > 4.0000e+05

cond(A, 2)
% > 4.0000e+05
cond(B, 2)
% > 4.0000e+05

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [2, 0, 5, 8; 0, 2, -1, -3; -2, 6, 2, -3; 4, -4, 0, 2];
[L U] = LU_nopvt(A);
L
% >  1.0000         0         0         0
% >  0         1.0000         0         0
% > -1.0000    3.0000    1.0000         0
% >  2.0000   -2.0000   -1.2000    1.0000
U
% >	 2.0000         0    5.0000    8.0000
% >	      0    2.0000   -1.0000   -3.0000
% >	      0         0   10.0000   14.0000
% >	      0         0         0   -3.2000


[L,U,P] = lu(A);

L
% >  1.0000         0         0         0
% > -0.5000    1.0000         0         0
% >  0.5000    0.5000    1.0000         0
% >  0         0.5000   -0.5000    1.0000

U
% > 4    -4     0     2
% > 0     4     2    -2
% > 0     0     4     8
% > 0     0     0     2

P
% > 0     0     0     1
% > 0     0     1     0
% > 1     0     0     0
% > 0     1     0     0


[L U] = lu(A);
b = [0, 1, 0 ,0]';
y = L \ b;
x = U \ y;
x
% >  0.5000
% >  0.7500
% > -1.0000
% >  0.5000
