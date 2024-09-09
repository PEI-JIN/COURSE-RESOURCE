function [output] = autocorrelation_DRFFT64(x)
% autocorrelation_DRFFT64 (Real Autocorrelation Program by Calling DRFFT64)
% autocorrelation_DRFFT64(x) computes the autocorrelation of x[n]=n*(-1)^n,
% for n = 1, 2, ..., 32.

N = 32;
x_64 = [x zeros(1,N)];      % Zero data extended padding

[X1,X2] = DRFFT64(x_64,x_64);     % Calling DRFFT64 once
Rx = IFFT64(X1.*conj(X2));        % Calling IFFT64 IFFT program once
Rx_time_reverse = Rx(1+mod(0:-1:1-N,N));    % Time Reverse (Modulo N)
output = [real(Rx_time_reverse(2:N)), real(Rx(1:N))];

end