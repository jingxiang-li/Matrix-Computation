clear; close all; clc;
cd /home/jingxiang-li/Documents/Courses/CSCI-5304/CSCI-5304-HW/HW4

X = [4, 8; 9, 3; 2, 6; 7, 4; 8, 3];
y = [4; 8; 1; 4; 6];

[Q R] = qr(X);
Q = Q(:, 1:2);
R = R(1:2, :);
%% Here we use inv() since the dimension of R is very small
theta = inv(R) * (Q' * y);
alpha_0 = theta(1);
beta_0 = theta(2);
%% Make prediction
X_new = [3, 8];
prediction = X_new * theta;

alpha_0
% > alpha_0 =  0.77036

beta_0
% > beta_0 =  0.0093010

prediction
% > prediction =  2.3855



clear; close all; clc
% Read Data
lsidata;
[n, m] = size(A);
k = 50;

% Now Calculate A_50 Using svd
[U, S, V] = svds(A, k);
A_50 = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';

% q is the query vector in the sparse form
q_tmp = [1897, 1, 1; n, 1, 0];
q = spconvert(q_tmp);

% result_origin
result_origin = q' * A;
% Only 112th row has non-zero output, suggesting that the query result is the 112th document
headlines(112)
% > 112 Smoking Alters Brain Chemical

% result_reduced
result_reduced = q' * A_50;
[s index] = sort(result_reduced, 'descend');
headlines(index(1:10))

% > 492 Ultrasound Reveals Thinking Brain
% > 240 Enzyme Linked To Brain Aneurysm
% > 176 Scans Reveal Brain Reacting to Cocaine
% > 112 Smoking Alters Brain Chemical
% > 423 Diet Drugs Affect Brain Cells
% > 177 Leptin: Possible Diabetes Treatment?
% > 448 Diet Drugs Affect Brain Cells
% > 131 Obesity Hormone Regulates Blood Sugar
% > 392 Brain Chemicals Mimic Marijuana
% > 251 Newborn Brain Link to Mental Ills