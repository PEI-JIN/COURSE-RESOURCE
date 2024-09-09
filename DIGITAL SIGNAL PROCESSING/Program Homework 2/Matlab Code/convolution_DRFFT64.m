function [output] = convolution_DRFFT64(x,h)
% convolution_DRFFT64 (A Linear Convolution Program by Calling DRFFT64)
% convolution_DRFFT64(x,h) computes the convolution of two 32 real-data,
% which are x[n] = [3, 6, 9, ..., 96] and h[n] = 32 - 2n,
% for n = 1, 2, ..., 32.

N = 32;
x_64 = [x zeros(1,N)];  % Zero data extended padding
h_64 = [h zeros(1,N)];  % Zero data extended padding

[X,H] = DRFFT64(x_64,h_64);     % Calling DRFFT64 once
IFFT64_ouput = IFFT64(X.*H);    % Calling IFFT64 IFFT program once
output = real(IFFT64_ouput(1:2*N-1));

end