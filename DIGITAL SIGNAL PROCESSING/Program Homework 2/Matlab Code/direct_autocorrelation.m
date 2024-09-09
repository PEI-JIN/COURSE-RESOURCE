function [output] = direct_autocorrelation(x)
% direct_autocorrelation (A Direct Real Autocorrelation Program)
% direct_autocorrelation(x) computes the autocorrelation of x[n]=n*(-1)^n,
% for n = 1, 2, ..., 32.

m = length(x);  % Length of x[n]
Rx = zeros(1,2*m-1);
for i = 1:m
    Rx(i)=sum(x(i:m).*x(1:m-i+1));
end
Rx_time_reverse = Rx(1+mod(0:-1:1-m,m));    % Time Reverse (Modulo 32)
output = [Rx_time_reverse(2:m), Rx(1:m)];

end