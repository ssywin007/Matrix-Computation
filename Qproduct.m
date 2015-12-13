function p = Qproduct(u, x, t)
%% if t == 1, compute p = Q'x = H_n - 1 * ... * H_1 * x 
%% if t == 0, compute p = Qx = H_1 * ... * H_n - 1 * x 
 %%
n = size(u, 2) + 1;
p = x;
if t == 1
    for i = 1: n - 1
        tempux = u(:, i);
        p = p - 2 * tempux * (tempux' * p);
    end
end
if t == 0
    for i = 1: n - 1
        tempuy = u(:, n - i);
        p = p - 2 * tempuy * (tempuy' * p);
    end
end
