%% QR Algorithm to solve the eigenvalue eigenvector problem
function [eigVec eigVal iterNum resArray] = QRAlgo(A, iterMax, eps0)
	U = eye(size(A));
	iterNum = 1;
	resArray = zeros(iterMax, 1);
	A0 = A;
	[Q R] = qr(A);
	A = R * Q;
	U = U * Q;

	resArray(iterNum) = norm(diag(A0 * U - diag(diag(A)) * U), 1);
	while (resArray(iterNum) > eps0 && iterNum < iterMax)
		[Q R] = qr(A);
		A = R * Q;
		U = U * Q;		
		iterNum = iterNum + 1;
		resArray(iterNum) = norm(diag(A0 * U - diag(diag(A)) * U), 1);
	end
	
	eigVal = A;
	eigVec = U;
end