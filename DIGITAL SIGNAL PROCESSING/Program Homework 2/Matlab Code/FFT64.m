function [X] = FFT64(x)
% FFT64 (Decimation-in-time FFT Algorithm)
% FFT64(x) computes 64-point FFT of x and returns the result,
% where x is a complex 64-point vector.

N = 64;

% Use FFT32
G = FFT32(x(1:2:N));    % 32-point FFT (Even part of x)
H = FFT32(x(2:2:N));    % 32-point FFT (Odd part of x)

W_N = exp(-1i*2*pi/N.*(0:N/2-1));  % Complex sinusoidal components

% Butterfly
X(1:N/2) = G + H.*W_N;
X(N/2+1:N) = G - H.*W_N;

end