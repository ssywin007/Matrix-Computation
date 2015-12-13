function [Q, R] = QRclassicdecomposition(A)
%% classic method (Schimidt Orthogonalization) to compute the QR
%% decomposition of A, m by n matrix and n <= m with full column rank
%%
m = size(A, 1);
n = size(A, 2);
r = rank(A);
if r < n
    fprintf('The input matrix does not have full column rank!\n')
end
Im = eye(m);
Q = zeros(m, m);
R = zeros(m, n);
%compute Q1, Q = [Q1, Q2]
for i = 1: n
    R(1: i - 1, i) = Q(:, 1: i - 1)' * A(:, i);
    temp = A(:, i) - Q(:, 1: i - 1) * R(1: i - 1, i);
    R(i, i) = sqrt(temp' * temp);
    Q(:, i) = temp / R(i, i);
end
%guarantee R is upper triangular, actually, it should be
R = triu(R);
%compute Q2
for i = (n + 1): m
    temp = rand(m, 1);
    temp = temp / sqrt(temp' * temp);
    temp = (temp - Q(:, 1: i - 1) * Q(:, 1: i - 1)' * temp);
    Q(:, i) = temp / sqrt(temp' * temp);
end
%guarantee Q is orthogonal
if norm(Q' * Q - Im) > 1e-10
    [QQ, QR] = QRclassicdecomposition(Q);
    Q = QQ;
    R = QR * R;
    R = triu(R);
    if norm(Q' * Q - Im) > 1e-10
        fprintf('ERROR, the Q of classic QR decomposition is not unitary! ||QT * Q - Im|| = %f\n', norm(Q' * Q - Im))
    end
end
if norm(Q * R - A) > 1e-10
    fprintf('ERROR, the classic QR decomposition fails! ||Q * R - A|| = %f\n', norm(Q * R - A))
end
