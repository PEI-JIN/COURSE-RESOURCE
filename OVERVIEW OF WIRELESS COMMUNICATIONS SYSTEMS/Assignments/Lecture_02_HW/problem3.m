% Problem # 03 (Matlab Assignment)
% Perform simulation and plot the bit error rate (BER) vs Eb/N0 (dB)
% of a BPSK communication system over Rayleigh fading and AWGN channel.
% Written by P.-J. Su 2022/11/08
% Run: To compile and run this computer program using MATLAB.

% Clear Command Window
clc
% Remove all variables from the current workspace,
% releasing them from system memory.
clear
% Close all figures whose handles are visible.
close all
% Start stopwatch timer.
tic

num_bits = 10^7;    % Number of data bit
Eb = 1;     % Energy per bit
EbOverN0_dB = 0:1:35;     % Eb/N0 in dB
EbOverN0 = 10.^(EbOverN0_dB/10);     % Eb/N0
num_of_bit_error_AWGN = zeros(size(EbOverN0_dB));    % Number of bit error (AWGN)
num_of_bit_error_Rayleigh = zeros(size(EbOverN0_dB));    % Number of bit error (Rayleigh)
BER_theo_AWGN = zeros(size(EbOverN0_dB));
BER_theo_Rayleigh = zeros(size(EbOverN0_dB));

% Mapping to the signal constellation
BPSK = [-sqrt(Eb), sqrt(Eb)];

for Idx = 1: length(EbOverN0_dB)
    tx_bit = ceil(2.*rand(1, num_bits))-1;  %% Transmitted data bit
    tx_sym = BPSK(tx_bit(1:num_bits)+1);  %% Transmitted symbol
    
    % AWGN channel
    sigma = sqrt(Eb/(2*EbOverN0(Idx)));
    n = randn(1, length(tx_sym))*sigma;
    
    % Rayleigh Fading channel
    x = randn(1, length(tx_sym))*sqrt(0.5);
    y = randn(1, length(tx_sym))*sqrt(0.5);
    h = sqrt(x.*x + y.*y);

    % The received signal at the detector (AWGN channel)
    rx_sym_AWGN = tx_sym + n;  
    
    rx_bit_AWGN = zeros(size(tx_bit)); %% Reset rx bit decisions 
    rx_bit_AWGN(rx_sym_AWGN > 0) = 1;   %% Decision of bit '1'
    err_pat_AWGN = xor(tx_bit, rx_bit_AWGN); %% Compare tx and rx bits

    % The received signal at the detector (Rayleigh)
    rx_sym_Rayleigh = h.*tx_sym + n;  
    
    rx_bit_Rayleigh = zeros(size(tx_bit)); %% Reset rx bit decisions 
    rx_bit_Rayleigh(rx_sym_Rayleigh > 0) = 1;   %% Decision of bit '1'
    err_pat_Rayleigh = xor(tx_bit, rx_bit_Rayleigh); %% Compare tx and rx bits
    
    % Error counting
    num_of_bit_error_AWGN(Idx) = num_of_bit_error_AWGN(Idx) + sum(err_pat_AWGN);
    num_of_bit_error_Rayleigh(Idx) = num_of_bit_error_Rayleigh(Idx) + sum(err_pat_Rayleigh);

    % Error probability calculation (Analysis)
    BER_theo_AWGN(Idx) = qfunc(sqrt(2*EbOverN0(Idx)));
    SNR = 2*EbOverN0(Idx);
    BER_theo_Rayleigh(Idx) = 0.5*(1-sqrt(SNR/(2+SNR)));
end

% Error probability calculation (Simulation)
BER_sim_AWGN = num_of_bit_error_AWGN/num_bits;
BER_sim_Rayleigh = num_of_bit_error_Rayleigh/num_bits;

% Plot the bit error rate (BER) curves
semilogy(EbOverN0_dB, BER_sim_AWGN, '-r');
hold on
semilogy(EbOverN0_dB, BER_theo_AWGN,'bo');
hold on
semilogy(EbOverN0_dB, BER_sim_Rayleigh, '-m');
hold on
semilogy(EbOverN0_dB, BER_theo_Rayleigh,'gd');
axis([0 35 1e-5 1]);
title('BER performance of a BPSK system');
legend('AWGN channel (Simulation)','AWGN channel (Analysis)', ...
       'Rayleigh Fading and AWGN channel (Simulation)', ...
       'Rayleigh Fading and AWGN channel (Analysis)');
xlabel('\itE_{b}\rm/\itN\rm_{0} (dB)');
ylabel('BER');
grid on
