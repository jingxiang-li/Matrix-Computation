close all; clear; clc;
cd /home/jingxiang-li/Documents/Courses/CSCI-5304/CSCI-5304-HW/HW3

m = 3;
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

plot(m_array, result_time)
hold on;
plot(m_array, result_time, 'o')
xlabel('m');ylabel('elapsed time');title('elapsed time vs. m');
fontsize=20;
set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
set([gca; findall(gca, 'Type','line')], 'linewidth', 3);
saveas(1, 'presentation.pdf');


%% Compare
m_array = 500 : 250 : 2500;
result_time = zeros(size(m_array)(2), 2);
result_res = zeros(size(m_array)(2), 2);
for i = 1 : size(m_array)(2)
    m = m_array(i);
    A_trid = [ [nan; ones(2*m,1)*m ], (m+1+(-m:m)') , [ones(2*m,1)*m;nan]];
    A_full = diag(A_trid(:,2))+diag(A_trid(2:end,1),-1)+diag(A_trid(1:end-1,3),1);
    b = ones(size(A_trid, 1), 1);
    tic();
        x1 = GE_Band(A_trid, b);
    result_time(i, 1) = toc();
    result_res(i, 1) = norm(b - A_full * x1, 1);
    tic();
    	x2 = A_full \ b;
    result_time(i, 2) = toc();
    result_res(i, 2) = norm(b - A_full * x2, 1);
endfor

[result_time result_res]
% > 7.6898e-02   1.0637e-01   1.6211e-12   1.9318e-14
% > 1.0888e-01   2.9426e-01   3.7359e-13   3.5971e-14
% > 1.5138e-01   7.2223e-01   7.7149e-13   5.1292e-14
% > 1.8058e-01   1.2491e+00   1.2113e-12   1.6820e-13
% > 2.1552e-01   2.0951e+00   3.8237e-12   5.5300e-13
% > 2.5187e-01   3.2953e+00   2.3949e-12   1.0525e-13
% > 2.8727e-01   4.8834e+00   1.0439e-12   8.1157e-14
% > 3.2394e-01   6.8370e+00   1.0638e-12   8.9595e-14
% > 3.7854e-01   9.3374e+00   2.6795e-12   1.2967e-13


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc

m = 3;
A_trid = [ [nan; ones(2*m,1)*m ], (m+1+(-m:m)') , [ones(2*m,1)*m;nan]];
A_full = diag(A_trid(:,2))+diag(A_trid(2:end,1),-1)+diag(A_trid(1:end-1,3),1);
b = ones(size(A_trid, 1), 1);

[R b1 Q] = QR_Given(A_trid, b);
x = UpTri_Band4(R, b1);
[A_full \ b x]
Q

n = size(A_trid, 1);
Q1 = eye(n);
for k = 1 : size(Q, 1)
    c = Q(k, 1);
    s = Q(k, 2);
    i = Q(k, 3);
    j = Q(k, 4);
    tmp = eye(n);
    tmp(i, i) = c;
    tmp(j, j) = c;
    tmp(i, j) = s;
    tmp(j, i) = -s;
    Q1 = Q1 * tmp;
end

m = 2;
A_trid = [ [nan; ones(2*m,1)*m ], (m+1+(-m:m)') , [ones(2*m,1)*m;nan]];
A_full = diag(A_trid(:,2))+diag(A_trid(2:end,1),-1)+diag(A_trid(1:end-1,3),1);
b = ones(size(A_trid, 1), 1);
n = size(A_trid, 1);
[R b1 Q] = QR_Given(A_trid, b);
Q = Q_trans(Q, n);
[Q1 R1] = qr(A_full);
% n = size(A_trid, 1);
% Q1 = eye(n);
% for k = 1 : size(Q, 1)
%     c = Q(k, 1);
%     s = Q(k, 2);
%     i = Q(k, 3);
%     j = Q(k, 4);
%     tmp = eye(n);
%     tmp(i, i) = c;
%     tmp(j, j) = c;
%     tmp(i, j) = s;
%     tmp(j, i) = -s;
%     Q1 = Q1 * tmp;
% end







Q
% > 0.44721   0.36515  -0.58321   0.43004   0.37629
% > 0.89443  -0.18257   0.29161  -0.21502  -0.18814
% > 0.00000   0.91287   0.29161  -0.21502  -0.18814
% > 0.00000   0.00000   0.69985   0.53755   0.47036
% > 0.00000   0.00000   0.00000   0.65850  -0.75258

Q1
% > -0.44721   0.36515   0.58321  -0.43004  -0.37629
% > -0.89443  -0.18257  -0.29161   0.21502   0.18814
% > -0.00000   0.91287  -0.29161   0.21502   0.18814
% > -0.00000   0.00000  -0.69985  -0.53755  -0.47036
% > -0.00000   0.00000  -0.00000  -0.65850   0.75258







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
