clear all;
close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
iteration = 400;
num_symbol = 200;

SNR_dB = 0:2:20; % in dB
SNR = 10.^(SNR_dB./10);

errors_ML = zeros(length(SNR_dB),1);
qam_symbol =0:1:sym_QAM-1;
qam_signal = qammod(qam_symbol, sym_QAM, 'UnitAveragePower', true); 
%%ML
%% s_ML = arg max p(s|y)
txsymbol_QAM = randi([0,sym_QAM-1], Ntx, num_symbol); % 生成0到15之间的整数作为符号
txsignal_QAM = qammod(txsymbol_QAM, sym_QAM, 'UnitAveragePower', true); % 使用灰度映射调制符号
H = sqrt(0.5)*(randn(Nrx,Ntx)+1i*randn(Nrx,Ntx));
for l = 1:length(SNR_dB)
    N0 = 1/(10^(SNR_dB(l)/10));
    for iter = 1:iteration
        
        noise = sqrt(N0/2)*(randn(Nrx,num_symbol)+1i*randn(Nrx,num_symbol));
        
        y = H * txsignal_QAM + noise;

        x_hat = zeros(Ntx, num_symbol);

        for index_symbol = 1: num_symbol
            Min = inf;
            for i = 1:sym_QAM
                for j = 1 : sym_QAM
                    temp_x_signal = [qam_signal(i);qam_signal(j)];
                    temp = norm(y(:,index_symbol)-H*temp_x_signal);
                    if temp < Min
                        x_hat(:,index_symbol) = temp_x_signal;
                        Min = temp;
                    end

                end

            end

        end
        QAM_demod_ML = qamdemod(x_hat,sym_QAM, 'UnitAveragePower', true);
        errors_ML(l) = errors_ML(l) + sum(sum(QAM_demod_ML~=txsymbol_QAM));
        
    end
    disp(l)
end
error_rate_ML = errors_ML/(iteration*num_symbol*Ntx);

figure;
semilogy(SNR_dB,error_rate_ML,'r-o');
legend('ML');
title('BER of ML');
xlabel('SNR in dB');
ylabel('BER');