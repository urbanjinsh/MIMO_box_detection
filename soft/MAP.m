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
QAM_demod_soft = zeros(Ntx,num_symbol);
QAM_demod_MMSE = zeros(Ntx,num_symbol);

errors_soft = zeros(length(SNR_dB),1);
errors_MMSE = zeros(length(SNR_dB),1);
noise_variance = zeros(length(SNR_dB),1);
txsymbol_QAM = randi([0,sym_QAM-1], Ntx, num_symbol); % 生成0到15之间的整数作为符号
txsignal_QAM = qammod(txsymbol_QAM, sym_QAM, 'UnitAveragePower', true); % 使用灰度映射调制符号
H = sqrt(0.5)*(randn(Nrx,Ntx)+1i*randn(Nrx,Ntx));

for l = 1:length(SNR_dB)
    N0 = 1/(10^(SNR_dB(l)/10));
    for iter = 1:iteration
        %%生成信号+调制

        w_MMSE = (H'*H+N0*eye(Ntx))\H';


        noise = sqrt(N0/2)*(randn(Nrx,num_symbol)+1i*randn(Nrx,num_symbol));
        y = H * txsignal_QAM + noise; 

        %hard
        y_MMSE = w_MMSE * y;

        %soft
        
        soft = soft_decision_2x2(y_MMSE);
        % SINR1 = (abs((w_MMSE(1,:)*H(:,1)))^2)/(abs((w_MMSE(1,:)*H(:,2)))^2+w_MMSE(1,:)*w_MMSE(1,:)'*N0);
        % SINR2 = (abs((w_MMSE(2,:)*H(:,2)))^2)/(abs((w_MMSE(2,:)*H(:,1)))^2+w_MMSE(2,:)*w_MMSE(2,:)'*N0);
        % QAM_demod_soft = soft .* [SINR1;SINR2];
        QAM_demod_soft = soft;
        



        %解调
       
        QAM_demod_MMSE = qamdemod(y_MMSE,sym_QAM, 'UnitAveragePower', true);

        %计算误码
        errors_soft(l) = errors_soft(l) + sum(sum(QAM_demod_soft~=txsymbol_QAM));
        errors_MMSE(l) = errors_MMSE(l) + sum(sum(QAM_demod_MMSE~=txsymbol_QAM));
    end

end

%计算误码率
error_rate_soft = errors_soft/(iteration*num_symbol*Ntx);
error_rate_MMSE = errors_MMSE/(iteration*num_symbol*Ntx);

figure;
semilogy(SNR_dB,error_rate_soft,'r-o',SNR_dB,error_rate_MMSE,'b-s');
legend('soft','MMSE');
title('BER of soft and MMSE');
xlabel('SNR in dB');
ylabel('BER');

function [x4_soft] = soft_decision_2x2(x)
x=x(:).'; 
X_r=real(x); 
X_i=imag(x);

d = 1/sqrt(10);
if X_i < -2*d
    X_1 = 2*d + 2*X_i;
elseif X_i > 2*d
    X_1 = -2*d + 2*X_i;
else
    X_1 = -X_i;
end

if X_r < -2*d
    X_3 = 2*d + 2*X_r;
elseif X_r > 2*d
    X_3 = -2*d + 2*X_r;
else 
    X_3 = X_r;
end
X_2 = 2*d-abs(X_r);
X_0 = 2*d-abs(X_i);
X=[X_3; X_2;  X_1; X_0];  
X_int_bin = X >0;
bin_to_dec = [2^3;2^2;2^1;2^0];
X_int_dec = X_int_bin.' * bin_to_dec;
x4_soft = reshape((X_int_dec(:)),2,[]);
end

