function [output] = direct_convolution(x,h)
% direct_convolution (A Direct Linear Convolution Program)
% direct_convolution(x,h) computes the convolution of two 32 real-data,
% which are x[n] = [3, 6, 9, ..., 96] and h[n] = 32 - 2n,
% for n = 1, 2, ..., 32.

m = length(x);  % Length of x[n]
n = length(h);  % Length of h[n]
output = zeros(1,m+n-1);
for i = 1:m+n-1
    for j = 1:m
        if(i-j+1 >= 1 && i-j+1 <= n)
             output(i) = output(i) + x(j)*h(i-j+1);
        end
    end
end

end