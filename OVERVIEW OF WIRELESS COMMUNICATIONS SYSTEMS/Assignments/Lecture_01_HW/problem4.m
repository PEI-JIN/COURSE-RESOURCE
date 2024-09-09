% Problem # 04 (Matlab Assignment)
% Perform simulation and plot the bit error rate (BER) vs Eb/N0 (dB)
% and symbol error rate (SER) vs Eb/N0 (dB) of a M-ary QAM communication
% system (M = 4, 16, 64) over AWGN channel.
% Written by P.-J. Su 2022/9/27
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

SNR_in_dB = 0:1:20;     % SNR (Eb/N0) in dB

 % 4-QAM
[BER_sim_1, BER_theo_1, SER_sim_1, SER_theo_1] = QAM(4, SNR_in_dB);
 % 16-QAM
[BER_sim_2, BER_theo_2, SER_sim_2, SER_theo_2] = QAM(16, SNR_in_dB);
 % 64-QAM
[BER_sim_3, BER_theo_3, SER_sim_3, SER_theo_3] = QAM(64, SNR_in_dB);

% Plot the bit error rate (BER) curves
clf
semilogy(SNR_in_dB, BER_sim_1, '-m');
hold on
semilogy(SNR_in_dB, BER_theo_1,'o','Color',[0 0.447 0.741]);
hold on
semilogy(SNR_in_dB, BER_sim_2, '--r');
hold on
semilogy(SNR_in_dB, BER_theo_2,'*','Color',[0.9290 0.6940 0.1250]);
hold on
semilogy(SNR_in_dB, BER_sim_3, '-g');
hold on
semilogy(SNR_in_dB, BER_theo_3,'^','Color',[0.4940 0.1840 0.5560]);
axis([0 20 1e-5 1]);
title('Bit error rate of a M-ary QAM system');
legend('4-QAM (Simulation)','4-QAM (Analysis)', ...
       '16-QAM (Simulation)','16-QAM (Analysis)', ...
       '64-QAM (Simulation)','64-QAM (Analysis)');
xlabel('\itE_{b}\rm/\itN\rm_{0} (dB)');
ylabel('BER');
grid on

% Plot the symbol error rate (SER) curves
figure
semilogy(SNR_in_dB, SER_sim_1, '-m');
hold on
semilogy(SNR_in_dB, SER_theo_1,'o','Color',[0 0.447 0.741]);
hold on
semilogy(SNR_in_dB, SER_sim_2, '--r');
hold on
semilogy(SNR_in_dB, SER_theo_2,'*','Color',[0.9290 0.6940 0.1250]);
hold on
semilogy(SNR_in_dB, SER_sim_3, '-g');
hold on
semilogy(SNR_in_dB, SER_theo_3,'^','Color',[0.4940 0.1840 0.5560]);
axis([0 20 1e-5 1]);
title('Symbol error rate of a M-ary QAM system');
legend('4-QAM (Simulation)','4-QAM (Analysis)', ...
       '16-QAM (Simulation)','16-QAM (Analysis)', ...
       '64-QAM (Simulation)','64-QAM (Analysis)');
xlabel('\itE_{b}\rm/\itN\rm_{0} (dB)');
ylabel('SER');
grid on

% Read elapsed time from stopwatch.
toc

function [BER_sim, BER_theo, SER_sim, SER_theo] = QAM(M, SNR_in_dB)
    N = 10^6;    % Number of data symbol
    k = log2(M);    % M = 2^k
    SNR = 10.^(SNR_in_dB/10);     % SNR (Eb/N0)
    dmin = 1;      % Minimum distance between symbol
    Es = 2*(M-1)*(dmin^2)/3;     % Energy per symbol
    num_of_symbol_error = zeros(size(SNR_in_dB));    % Number of symbol error
    num_of_bit_error = zeros(size(SNR_in_dB));    % Number of bit error
    SER_theo = zeros(size(SNR_in_dB));
    BER_theo = zeros(size(SNR_in_dB));
    
    % Generation of transmitted bits based on the data source (Gray code)
    QAM_4_Gray_code = [0 1; 1 1; 0 0; 1 0];
    QAM_16_Gray_code = [0 0 1 0; 0 1 1 0; 1 1 1 0; 1 0 1 0; ...
                        0 0 1 1; 0 1 1 1; 1 1 1 1; 1 0 1 1; ...
                        0 0 0 1; 0 1 0 1; 1 1 0 1; 1 0 0 1; ...
                        0 0 0 0; 0 1 0 0; 1 1 0 0; 1 0 0 0];
    QAM_64_Gray_code = [0 0 0 1 0 0; 0 0 1 1 0 0; 0 1 1 1 0 0; 0 1 0 1 0 0; ...
                        1 1 0 1 0 0; 1 1 1 1 0 0; 1 0 1 1 0 0; 1 0 0 1 0 0; ...
                        0 0 0 1 0 1; 0 0 1 1 0 1; 0 1 1 1 0 1; 0 1 0 1 0 1; ...
                        1 1 0 1 0 1; 1 1 1 1 0 1; 1 0 1 1 0 1; 1 0 0 1 0 1; ...
                        0 0 0 1 1 1; 0 0 1 1 1 1; 0 1 1 1 1 1; 0 1 0 1 1 1; ...
                        1 1 0 1 1 1; 1 1 1 1 1 1; 1 0 1 1 1 1; 1 0 0 1 1 1; ...
                        0 0 0 1 1 0; 0 0 1 1 1 0; 0 1 1 1 1 0; 0 1 0 1 1 0; ...
                        1 1 0 1 1 0; 1 1 1 1 1 0; 1 0 1 1 1 0; 1 0 0 1 1 0; ...
                        0 0 0 0 1 0; 0 0 1 0 1 0; 0 1 1 0 1 0; 0 1 0 0 1 0; ...
                        1 1 0 0 1 0; 1 1 1 0 1 0; 1 0 1 0 1 0; 1 0 0 0 1 0; ...
                        0 0 0 0 1 1; 0 0 1 0 1 1; 0 1 1 0 1 1; 0 1 0 0 1 1; ...
                        1 1 0 0 1 1; 1 1 1 0 1 1; 1 0 1 0 1 1; 1 0 0 0 1 1; ...
                        0 0 0 0 0 1; 0 0 1 0 0 1; 0 1 1 0 0 1; 0 1 0 0 0 1; ...
                        1 1 0 0 0 1; 1 1 1 0 0 1; 1 0 1 0 0 1; 1 0 0 0 0 1; ...
                        0 0 0 0 0 0; 0 0 1 0 0 0; 0 1 1 0 0 0; 0 1 0 0 0 0; ...
                        1 1 0 0 0 0; 1 1 1 0 0 0; 1 0 1 0 0 0; 1 0 0 0 0 0];
    
    % Mapping to the signal constellation
    QAM_mapping = zeros(M,2);
    QAM_mapping_Idx = 1;
    for i = 1:sqrt(M)    
        for j = 1:sqrt(M)
            QAM_mapping(QAM_mapping_Idx,1) = -(sqrt(M)-1) + 2*(j-1);
            QAM_mapping(QAM_mapping_Idx,2) = (sqrt(M)-1) - 2*(i-1);
            QAM_mapping_Idx = QAM_mapping_Idx + 1;
        end
    end
    
    for snr_Idx = 1:length(SNR_in_dB)
        data_source = zeros(N,1);   % Data source
        tx_signal = zeros(N,2);     % Transmitted symbols
        rx_signal = zeros(N,2);     % Received symbols
        tx_bit = zeros(N*k,1);      % Transmitted bits
        rx_bit = zeros(N*k,1);      % Received bits
        
        for i = 1:N
            % Generation of the data source
            data_source(i) = floor(M*rand) + 1;
    
            % The transmitted signal at the M-ary PSK mapper
            tx_signal(i,:) = QAM_mapping(data_source(i),:);
        
            % AWGN channel
            sigma = sqrt(Es/(2*k*SNR(snr_Idx)));
            noise = randn(1,2)*sigma;
        
            % The received signal at the detector
            rx_signal(i,:) = tx_signal(i,:) + noise;
    
            % Metric computation
            metrics = zeros(M,1);
            % Detection and symbol error calculation
            for j = 1:M
                metrics(j) = (rx_signal(i,1)-QAM_mapping(j,1))^2 + ...
                             (rx_signal(i,2)-QAM_mapping(j,2))^2;
            end
            [~, decision] = min(metrics); 
            if (decision ~= data_source(i))
                num_of_symbol_error(snr_Idx) = num_of_symbol_error(snr_Idx)+1;
            end
           
             % Generation of transmitted bits and received bits
            for l = 1:k
                switch M
                    case 4 
                        tx_bit((i-1)*k + l) = QAM_4_Gray_code(data_source(i),l);
                        rx_bit((i-1)*k + l) = QAM_4_Gray_code(decision,l);
                    case 16 
                        tx_bit((i-1)*k + l) = QAM_16_Gray_code(data_source(i),l);
                        rx_bit((i-1)*k + l) = QAM_16_Gray_code(decision,l);
                    case 64 
                        tx_bit((i-1)*k + l) = QAM_64_Gray_code(data_source(i),l);
                        rx_bit((i-1)*k + l) = QAM_64_Gray_code(decision,l);
                end    
            end
        end
    
        % Bit error calculation (Simulation)
        err_pat = xor(tx_bit, rx_bit);      % compare tx and rx bits
        num_of_bit_error(snr_Idx) = num_of_bit_error(snr_Idx) + sum(err_pat);
    
        % Error probability calculation (Analysis)
        Q = qfunc(dmin/sigma);
        const = 1/sqrt(M) - 1;
        SER_theo(snr_Idx) = -4*const*Q*(1 + const*Q);    % Theoretical SER
        BER_theo(snr_Idx) = SER_theo(snr_Idx)/k;    % Theoretical BER
    end
    
    % Error probability calculation (Simulation)
    SER_sim = num_of_symbol_error/(N);   % Simulated SER
    BER_sim = num_of_bit_error/(N*k);    % Simulated BER
end