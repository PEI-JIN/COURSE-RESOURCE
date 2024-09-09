function [X,Y] = DRFFT64(x,y)
% DRFFT64 (Double-real 64-point FFT Program)
% DRFFT64(x) is a double-real 64 FFT program,
% where x and y are two real 64-point data vectors.

N = 64;

z = x + y*1i;   % z[n] = x[n] + j*y[n]

Z = FFT64(z);     % Calling FFT64 FFT program once

X = 0.5*(Z+conj(Z(1+mod(0:-1:1-N,N))));     % Symmetry Property: Re{z[n]}

Y = 0.5*(Z-conj(Z(1+mod(0:-1:1-N,N))))/1i;  % Symmetry Property: j*Im{z[n]}

end