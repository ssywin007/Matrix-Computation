function u = u_householder(x, m)
%% calculate the u of householder decomposition, Hu * x = (I - 2 * u * u') *
%%  x= +- || x || * e1, m is the dimension of u, if dim(x) < m, set first m
%%  -dim(x) elements of u 0
%%
u = zeros(m, 1);
n = size(x, 1);
% set the sigh of alpha
if x(1) < 0
    alpha = sqrt(x' * x);
else
    alpha = - sqrt(x' * x);
end
u1 = sqrt((1 - x(1) / alpha) / 2);
u(m - n + 1) = u1;
if n > 1
    u((m - n + 2) : m) = -1 * x(2 : n) / (2 * alpha * u1);
end
% u should be unit vectors
while (u' * u) - 1 > 1e-10
    u((m - n + 1) : m) = u((m - n + 1) : m) / sqrt(u' * u);
end

