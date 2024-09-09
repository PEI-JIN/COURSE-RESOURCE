%%% DSP Homework #1, Written by SU PEI-JIN 2022/12/27
% Clear Command Window
clc
% Remove all variables from the current workspace,
% releasing them from system memory.
clear
% Close all figures whose handles are visible.
close all

% Generator initialization
seed = 1;
% Initialize the generator using a seed of 1.
rng(seed);

%% Verification of a complex 32-point FFT program
disp('Verification of a complex 32-point FFT program');
% x is a complex 32-point vector
N = 32;
x = randn(1,N) + 1i*randn(1,N);

% Start stopwatch timer.
tic

% Use MATLAB built-in function 
X_fft = fft(x);
% Use FFT32
X_FFT32 = FFT32(x);

% Create figure window.
figure('Name','Verification of a complex 32-point FFT program');
stem(abs(X_fft),'-or')
hold on;
stem(abs(X_FFT32),'--*b')
title('Verification of a complex 32-point FFT program: X[n]');
legend('fft (MATLAB built-in)', 'FFT32(x)');
% Label x-axis.
xlabel('n');

% Read elapsed time from stopwatch.
toc

%% Verification of a complex 64-point FFT program
disp('Verification of a complex 64-point FFT program');
% x is a complex 64-point vector
N = 64;
x = randn(1,N) + 1i*randn(1,N);

% Start stopwatch timer.
tic

% Use MATLAB built-in function 
X_fft = fft(x);
% Use FFT64
X_FFT64 = FFT64(x);

% Create figure window.
figure('Name','Verification of a complex 64-point FFT program');
stem(abs(X_fft),'-or')
hold on;
stem(abs(X_FFT64),'--*b')
title('Verification of a complex 64-point FFT program: X[n]');
legend('fft (MATLAB built-in)', 'FFT64(x)');
% Label x-axis.
xlabel('n');

% Read elapsed time from stopwatch.
toc

%% Verification of a double-real 64-point FFT program
disp('Verification of a double-real 64-point FFT program');
% x and y are two real 64-point data vector
N = 64;
x = randn(1,N);
y = randn(1,N);

% Start stopwatch timer.
tic

% Use MATLAB built-in function 
X_fft = fft(x);
Y_fft = fft(y);
% Use DRFFT64
[X_DRFFT64,Y_DRFFT64] = DRFFT64(x,y);

% Create figure window.
figure('Name','Verification of a double-real 64-point FFT program: X[n]');
stem(abs(X_fft),'-or')
hold on;
stem(abs(X_DRFFT64),'--*b')
title('Verification of a double-real 64-point FFT program: X[n]');
legend('fft (MATLAB built-in)', 'DRFFT64(x,y)');
% Label x-axis.
xlabel('n');

figure('Name','Verification of a double-real 64-point FFT program: Y[n]');
stem(abs(Y_fft),'-dr')
hold on;
stem(abs(Y_DRFFT64),'--+b')
title('Verification of a double-real 64-point FFT program: Y[n]');
legend('fft (MATLAB built-in)', 'DRFFT64(x,y)');
% Label x-axis.
xlabel('n');

% Read elapsed time from stopwatch.
toc

%% Verification of an inverse complex 64-point FFT program
disp('Verification of an inverse complex 64-point FFT program');
% Compute the 64-point FFT of x, where x is a complex 64-point vector
N = 64;
x = randn(1,N) + 1i*randn(1,N);
X = fft(x);

% Start stopwatch timer.
tic

% Use MATLAB built-in function 
x_ifft = ifft(X);
% Use IFFT64
x_IFFT64 = IFFT64(X);

% Create figure window.
figure('Name','Verification of an inverse complex 64-point FFT program');
stem(abs(x_ifft),'-or')
hold on;
stem(abs(x_IFFT64),'--*b')
title('Verification of an inverse complex 64-point FFT program');
legend('ifft (MATLAB built-in)', 'IFFT64(X)');
% Label x-axis.
xlabel('n');

% Read elapsed time from stopwatch.
toc
