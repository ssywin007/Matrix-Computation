function [uh, R] = QRhouseholderdecomposition_u(A)
%% householder method to compute the QR decomposition of A, m by n matrix and n <= m with full column rank
%% return uh and R, A = QR and Q = H1 * H2 * ... * Hn, Hi = I - 2 * uhi * uhi'
%%
m = size(A, 1);
n = size(A, 2);
r = rank(A);
if r < n
    fprintf('The input matrix does not have full column rank!\n')
end
R = A;
uh = zeros(m, n);
%store u and compute R
for i = 1: n
     tempu = u_householder(R(i : m, i), m);
     uh(:, i) = tempu;
     R = R - 2 * tempu * (tempu' * R);
end
R = triu(R);
