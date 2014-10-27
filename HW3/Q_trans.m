%% Q_trans Transform Q from [c, s, i, j] to matrix
%% Input    Q from QR_Given
%%          n size of Given transformation matrix
%% functionname: function description

function [Q1] = Q_trans(Q, n)

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

end