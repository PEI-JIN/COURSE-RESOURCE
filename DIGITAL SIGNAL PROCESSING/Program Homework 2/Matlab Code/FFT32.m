function [X] = FFT32(x)
% FFT32 (Decimation-in-frequency FFT Algorithm)
% FFT32(x) computes 32-point FFT of x and returns the result,
% where x is a complex 32-point vector.

N = 32;

W_N = exp(-1i*2*pi./2.^(log2(N):-1:1));  % Complex sinusoidal components

for m = 1:log2(N)   % m = log2(N) stages (N-point FFT)
    for k = 1:2^(m-1)   % Each stage with k = 2^(m-1) parts
        for btf_Num = 1:N/(2^m)     % Each part needs N/(2^m) butterflies
            btf_Idx_1 = (k-1)*N/(2^(m-1)) + btf_Num;            
            btf_Idx_2 = (k-1)*N/(2^(m-1)) + btf_Num + N/(2^m); 
            X1 = x(btf_Idx_1);
            X2 = x(btf_Idx_2);
            x(btf_Idx_1) = X1 + X2;
            x(btf_Idx_2) = (X1 - X2)*(W_N(m)^(btf_Num-1));
        end
    end
end

X = x(bitrevorder((1:N)-1)+1);  % Permute data into bit-reversed order

end