s = RandStream('mt19937ar','Seed', 1994);
RandStream.setGlobalStream(s);

%relative errors
errorRc = zeros(499, 1);
errorRh = zeros(499, 1);
errorQc = zeros(499, 1);
errorQh = zeros(499, 1);
errorAc = zeros(499, 1);
errorAh = zeros(499, 1);
for n = 2:500
    B = randn(n, n);
    %guarantee B has full rank
    while rank(B) < n
        B = randn(n, n);
    end
    R = randn(n, n);
    %guarantee R has full rank
    while all(diag(R)) == 0
        R = randn(n, n);
    end
    R = triu(R);
    [QB, RB] = qr(B);
    A = QB * R;
    [QAc, RAc] = QRclassicdecomposition(A);
    [QAh, RAh] = QRhouseholderdecomposition(A);
    errorQc(n - 1) = norm(QB - QAc, 'fro');
    errorQh(n - 1) = norm(QB - QAh, 'fro');
    errorAc(n - 1) = norm(A - QAc * RAc, 'fro') / norm(A, 'fro');
    errorAh(n - 1) = norm(A - QAh * RAh, 'fro') / norm(A, 'fro');
    R = (diag(diag(R) > 0) + diag(diag(R) < 0) * -1) * R;
    RAh = (diag(diag(RAh) > 0) + diag(diag(RAh) < 0) * -1) * RAh;
    errorRc(n - 1) = norm(R - RAc, 'fro') / norm(R, 'fro');
    errorRh(n - 1) = norm(R - RAh, 'fro') / norm(R, 'fro');
end

scatter(2: 500, errorAc, 'k'); lsline; gtext('classic');
hold on; scatter(2: 500, errorAh, '*'); lsline; xlabel('n'); ylabel('relative error of A = QR'), title('A = QR - n'); gtext('householder');
scatter(2: 500, errorRc, 'k'); lsline; gtext('classic');
hold on; scatter(2: 500, errorRh, '*'); lsline; xlabel('n'); ylabel('relative error of R'), title('R - n'); gtext('householder');
scatter(2: 500, errorQc, 'k'); lsline; gtext('classic');
hold on; scatter(2: 500, errorQh, '*'); lsline; xlabel('n'); ylabel('relative error of Q'), title('Q - n'); gtext('householder');

scatter(errorAc, errorRc, 'k'); lsline; xlabel('relative error of A = QR'); ylabel('relative error of R'), title('R - A = QR, classic');
scatter(errorAh, errorRh, 'k'); lsline; xlabel('relative error of A = QR'); ylabel('relative error of R'), title('R - A = QR, householder');
scatter(errorAc, errorQc, '*'); lsline; xlabel('relative error of A = QR'); ylabel('error of Q'), title('Q - A = QR, classic');
scatter(errorAh, errorQh, '*'); lsline; xlabel('relative error of A = QR'); ylabel('error of Q'), title('Q - A = QR, householder');
 
scatter(errorAh,  errorAc); lsline; xlabel('relative error of QR in Householder'); ylabel('relative error of QR in Classic'), title('reletive error of QR, Classic- Householder'); 
   
   % vandermonde matrix
   m = 100;
   n = 15;
   alpha = linspace(0,m - 1,m)' / (m - 1);
   A = zeros(m, n);
   for i = 1: n 
       A(:, i) = alpha .^ (i - 1);
   end
   b = exp(sin(4 * alpha));
   % A = QR
   %householder
   [uA1, RA1] = QRhouseholderdecomposition_u(A);
   c = Qproduct(uA, b, 1);
   c14_1 = c(n) / RA(n, n);
   c14_1 - 2006.787453080206
   %classic
   [QA2, RA2] = QRclassicdecomposition(A);
   c = QA2' * b;
   c14_1 = c(n) / RA2(n, n);
   c14_1 - 2006.787453080206
   % [A, b] = QR
   %classic
   [Qaug1, Raug1] = QRclassicdecomposition([A, b]);
   c14_1 = Raug(15, 16) / Raug(15, 15);
   c14_1 - 2006.787453080206
   %householder
   [uaug2, Raug2] = QRhouseholderdecomposition([A, b]);
   c14_1 = Raug1(15, 16) / Raug1(15, 15);
   c14_1 - 2006.787453080206
   %normal equation
   [Qne, Rne] = QRclassicdecomposition(A' * A);     %A' * A doesn't have full col rank
   c = Qne' * A' * b;
   c14_1 = c(n) / Rne(n, n); 
   [Qne, Rne] = QRhouseholderdecomposition(A' * A);     %A' * A doesn't have full col rank
   c = Qne' * A' * b;
   c14_1 = c(n) / Rne(n, n); 
   %solve normal equation directly
   c = (A' * A) \ (A' * b);
   c(15)
