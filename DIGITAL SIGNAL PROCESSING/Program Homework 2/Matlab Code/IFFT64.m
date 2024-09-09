function [x] = IFFT64(X)
% IFFT64 (Inverse complex 64-point FFT Program)
% IFFT64(X) is an inverse complex 64-point FFT program,
% where X is a complex 64-point vectors.

N = 64;

x_tilde = FFT64(X)/N;    % Calling FFT64 FFT program once

x = x_tilde(1+mod(0:-1:1-N,N));     % Time Reverse (Modulo N)

end