%%% DSP Homework #2, Written by SU PEI-JIN 2022/12/28
% Clear Command Window
clc
% Remove all variables from the current workspace,
% releasing them from system memory.
clear
% Close all figures whose handles are visible.
close all

%% Verification of a real linear convolution program
disp('Verification of a real linear convolution program');

% x[n] = [3, 6, 9, ..., 96]
N = 32;
x = 3*(1:N);
% h[n] = 32 - 2n for n = 1, 2, ..., 32.
h = 32-2*(1:N);

% Start stopwatch timer.
tic

% Use MATLAB built-in function 
x_conv = conv(x,h);
% Use a direct real convolution program 
[x_direct_convolution] = direct_convolution(x,h);
% Use a real convolution program by calling DRFFT64(x,y) once
[x_convolution_DRFFT64] = convolution_DRFFT64(x,h);

% Create figure window.
figure('Name','Verification of a real linear convolution program (1)');
stem(x_conv,'-or')
hold on;
stem(x_direct_convolution,'--*b')
title('Verification of a direct real convolution program');
legend('conv (MATLAB built-in)', 'direct convolution');
% Label x-axis.
xlabel('n');

% Create figure window.
figure('Name','Verification of a real linear convolution program (2)');
stem(x_conv,'-or')
hold on;
stem(x_convolution_DRFFT64,'--*b')
title('Verification of a real convolution program (DRFFT64)');
legend('conv (MATLAB built-in)', 'convolution (DRFFT64)');
% Label x-axis.
xlabel('n');

% Create figure window.
figure('Name','Verification of a real linear convolution program (3)');
stem(x_conv,'-^r')
hold on;
stem(x_direct_convolution,'--*b')
hold on;
stem(x_convolution_DRFFT64,'-.ok')
title('Verification of a real linear convolution program');
legend('conv (MATLAB built-in)', 'direct convolution', ...
       'convolution (DRFFT64)');
% Label x-axis.
xlabel('n');

% Read elapsed time from stopwatch.
toc

%% Verification of a real data autocorrelation program
disp('Verification of a real data autocorrelation program');
% x[n]=n*(-1)^n for n = 1, 2, ..., 32.
N = 32;
x = (1:N).*(-1).^(1:N);

% Start stopwatch timer.
tic

% Use MATLAB built-in function 
x_xcorr = xcorr(x); % xcorr(x) returns the autocorrelation sequence of x.
% Use a direct real autocorrelation program
x_direct_autocorrelation = direct_autocorrelation(x);
% Use a real autocorrelation program by calling DRFFT64(x,y) once
x_autocorrelation_DRFFT64 = autocorrelation_DRFFT64(x);

% Create figure window.
figure('Name','Verification of a real data autocorrelation program (1)');
stem(-(N-1):(N-1),x_xcorr,'-or')
hold on;
stem(-(N-1):(N-1),x_direct_autocorrelation,'--*b')
title('Verification of a direct real autocorrelation program');
legend('xcorr (MATLAB built-in)', 'direct autocorrelation');
% Label x-axis.
xlabel('n');

% Create figure window.
figure('Name','Verification of a real data autocorrelation program (2)');
stem(-(N-1):(N-1),x_xcorr,'-or')
hold on;
stem(-(N-1):(N-1),x_autocorrelation_DRFFT64,'--*b')
title('Verification of a computation program (DRFFT64)');
legend('xcorr (MATLAB built-in)', 'autocorrelation (DRFFT64)');
% Label x-axis.
xlabel('n');

% Create figure window.
figure('Name','Verification of a real data autocorrelation program (3)');
stem(-(N-1):(N-1),x_xcorr,'-^r')
hold on;
stem(-(N-1):(N-1),x_direct_autocorrelation,'--*b')
hold on;
stem(-(N-1):(N-1),x_autocorrelation_DRFFT64,'-.ok')
title('Verification of a real data autocorrelation program');
legend('xcorr (MATLAB built-in)', 'direct autocorrelation', ...
       'autocorrelation (DRFFT64)');
% Label x-axis.
xlabel('n');

% Read elapsed time from stopwatch.
toc
