clear all;
close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
iteration = 400;
num_symbol = 200;
%EbN0s = 0:1:20;
SNR_dB = 0:2:20; % in dB
SNR = 10.^(SNR_dB./10);


%y_QAM = zeros(1,iteration);
QAM_demod_ZF = zeros(Ntx,num_symbol);
QAM_demod_MMSE = zeros(Ntx,num_symbol);

errors_ZF = zeros(length(SNR_dB),1);
errors_MMSE = zeros(length(SNR_dB),1);
noise_variance = zeros(length(SNR_dB),1);
txsymbol_QAM = randi([0,sym_QAM-1], Ntx, num_symbol); % 生成0到15之间的整数作为符号
txsignal_QAM = qammod(txsymbol_QAM, sym_QAM, 'UnitAveragePower', true); % 使用灰度映射调制符号
H = sqrt(0.5)*(randn(Nrx,Ntx)+1i*randn(Nrx,Ntx));

for l = 1:length(SNR_dB)
    N0 = 1/(10^(SNR_dB(l)/10));
    %omiga2 = 0.5*N0;
    for iter = 1:iteration
        %%生成信号+调制

        w_MMSE = (H'*H+N0*eye(Ntx))\H';
        %P_s = mean(abs(y(:)).^2);
        % y_noise = awgn(y,SNR_dB(l));

        % noise = y_noise - y;
        % noise_variance1 = var(noise(:));
        % sigPower = sum(abs(y(:)).^2)/length(y(:));
        % sigPower = mean(y.^2, 2);
        % noisePower = sigPower/SNR(l);
        % noise_variance = noisePower;

        noise = sqrt(N0/2)*(randn(Nrx,num_symbol)+1i*randn(Nrx,num_symbol));
        y = H * txsignal_QAM + noise;

        % w_ZF = inv(H'*H)*H';
        
        % w_MMSE = inv(H'*H+noise_variance*eye(Ntx))*H';
        
        y_MMSE = w_MMSE * y;

        %解调
        QAM_demod_ZF = qamdemod(y_ZF,sym_QAM, 'UnitAveragePower', true);
        QAM_demod_MMSE = qamdemod(y_MMSE,sym_QAM, 'UnitAveragePower', true);

        %计算误码
        errors_ZF(l) = errors_ZF(l) + sum(sum(QAM_demod_ZF~=txsymbol_QAM));
        errors_MMSE(l) = errors_MMSE(l) + sum(sum(QAM_demod_MMSE~=txsymbol_QAM));
    end

end

%计算误码率
error_rate_ZF = errors_ZF/(iteration*num_symbol*Ntx);
error_rate_MMSE = errors_MMSE/(iteration*num_symbol*Ntx);

figure;
semilogy(SNR_dB,error_rate_ZF,'r-o',SNR_dB,error_rate_MMSE,'b-s');
legend('ZF','MMSE');
title('BER of ZF and MMSE');
xlabel('SNR in dB');
ylabel('BER');

function [x4_soft] = soft_output2x2(x)

%MIMO-OFDM Wireless Communications with MATLAB㈢   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

sq10=sqrt(10); 
sq10_2=2/sq10;  
x=x(:).'; 
xr=real(x); 
xi=imag(x); 
X=[-xi; sq10_2-abs(xi); xr; sq10_2-abs(xr)]; 
x4_soft = X(:).'; 
end


