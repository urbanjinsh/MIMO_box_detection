clear all;
close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
k = log2(sym_QAM);
iteration = 100;
num_symbol = 100;
k_best_num = 4;

qam_symbol =0:1:sym_QAM-1;
qam_signal = qammod(qam_symbol, sym_QAM, 'UnitAveragePower', true); 

trellis = poly2trellis(7,[171 133]);
tbl = 32;
rate = 1/2;

%EbN0s = 0:1:20;
SNR_dB = 0:5:20; % in dB
SNR = 10.^(SNR_dB./10);


%y_QAM = zeros(1,iteration);
% QAM_demod_soft = zeros(Ntx,num_symbol);
QAM_demod_MMSE = zeros(Ntx,num_symbol);
kbest_decoder = zeros(Ntx,num_symbol);

errors_soft = zeros(length(SNR_dB),1);
errors_MMSE = zeros(length(SNR_dB),1);
errors_KB = zeros(length(SNR_dB),1);
noise_variance = zeros(length(SNR_dB),1);




for l = 1:length(SNR_dB)
    N0 = 1/(10^(SNR_dB(l)/10));
    for iter = 1:iteration
        %%生成信号+调制
        txsymbol_QAM = randi([0,sym_QAM-1],num_symbol,Ntx);
        txsymbol_QAM_bit = int_to_bin(txsymbol_QAM);
        
        txsymbol_QAM_enc1 = convenc(txsymbol_QAM_bit(:, 1),trellis);
        txsymbol_QAM_enc2 = convenc(txsymbol_QAM_bit(:, 2),trellis);
        txsymbol_QAM_enc = [txsymbol_QAM_enc1,txsymbol_QAM_enc2];
        txsignal_QAM = qammod(txsymbol_QAM_enc, sym_QAM,'InputType','bit', 'UnitAveragePower', true);
        % txsignal_QAM= txsignal_QAM';
        H = sqrt(0.5)*(randn(Nrx,Ntx)+1i*randn(Nrx,Ntx));

        w_MMSE = (H'*H+N0*eye(Ntx))\H';


        % noise = sqrt(N0/2)*(randn(Nrx,num_symbol)+1i*randn(Nrx,num_symbol));
        % y = H * txsignal_QAM' + noise; 
        y_withoutNoise = H * txsignal_QAM.';
        y = awgn(y_withoutNoise,SNR_dB(l));

        

        y_MMSE = w_MMSE * y;
        

        %解调
        x_hat_KB_soft = zeros(2*num_symbol*k, Ntx);
        index_symbol = 1;
        while index_symbol <= 2*num_symbol
            [best_path,r] = kbest_soft(H, y(:, index_symbol), qam_signal, k_best_num,N0);
            
            x_hat_KB_soft([k*index_symbol-k+1:k*index_symbol],: ) = r;
            index_symbol = index_symbol + 1;
        end
        % QAM_demod_soft = qamdemod(y_MMSE.',sym_QAM, 'OutputType','approxllr','UnitAveragePower', true, 'NoiseVariance', N0);
        % QAM_demod_soft = soft_decode(y_MMSE');
        QAM_demod_MMSE = qamdemod(y_MMSE.',sym_QAM,'OutputType','bit','UnitAveragePower', true);
        
        x_hat_KB = zeros(size(x_hat_KB_soft));
        x_hat_KB(x_hat_KB_soft<0) =1;

        symbol_hard_kbest1 = vitdec(x_hat_KB(:,1),trellis,tbl,'cont','hard');
        symbol_hard_kbest2 = vitdec(x_hat_KB(:,2),trellis,tbl,'cont','hard');

        symbol_hard1 = vitdec(QAM_demod_MMSE(:,1),trellis,tbl,'cont','hard');
        symbol_hard2 = vitdec(QAM_demod_MMSE(:,2),trellis,tbl,'cont','hard');
        symbol_soft1 = vitdec(x_hat_KB_soft(:,1),trellis,tbl,'cont','unquant');
        symbol_soft2 = vitdec(x_hat_KB_soft(:,2),trellis,tbl,'cont','unquant');

        symbol_hard = [symbol_hard1,symbol_hard2];
        symbol_hard_KB = [symbol_hard_kbest1,symbol_hard_kbest2];
        symbol_soft = [symbol_soft1,symbol_soft2];

        symbol_soft_dec = bin_to_dec(symbol_soft);
        symbol_hard_dec = bin_to_dec(symbol_hard);
        symbol_hard_KB_dec = bin_to_dec(symbol_hard_KB);

        symbol_soft_dec = symbol_soft_dec(tbl/k+1:end);
        symbol_hard_dec = symbol_hard_dec(tbl/k+1:end);
        symbol_hard_KB_dec = symbol_hard_KB_dec(tbl/k+1:end);

        %计算误码
        
        errors_soft(l) = errors_soft(l) + sum(sum(symbol_soft_dec~=txsymbol_QAM(1:end-tbl/k)));
        errors_MMSE(l) = errors_MMSE(l) + sum(sum(symbol_hard_dec~=txsymbol_QAM(1:end-tbl/k)));
        errors_KB(l) = errors_KB(l) + sum(sum(symbol_hard_KB_dec~=txsymbol_QAM(1:end-tbl/k)));
    end
    disp(l);

end

%计算误码率
error_rate_soft = errors_soft/(iteration*num_symbol*Ntx);
error_rate_MMSE = errors_MMSE/(iteration*num_symbol*Ntx);
error_rate_KB = errors_KB/(iteration*num_symbol*Ntx);

figure;
semilogy(SNR_dB,error_rate_soft,'r-o',SNR_dB,error_rate_MMSE,'b-s',SNR_dB,error_rate_KB,'g-s');
legend('soft','MMSE','KB');
title('BER of soft and MMSE');
xlabel('SNR in dB');
ylabel('BER');

function [binary_array] = int_to_bin(array)
[rows, cols] = size(array);

% 将数组中的每个元素转换为4位二进制字符串
binary_strings = dec2bin(array, 4);

% 创建一个2*4n的double型数组
binary_array = zeros(4 *rows,  cols);

% 将二进制字符串转换为数字，并存储在2*4n的double型数组中
for j = 1:cols
    for i = 1:rows
        % 获取当前元素的二进制字符串
        bin_str = binary_strings((j-1)*rows + i, :);
        % 将二进制字符串的每一位转换为数字，并存储在结果数组中
        for k = 1:4
            binary_array(4*(i-1) + k, j) = str2double(bin_str(k));
        end
    end
end
end

function [decimal_array] = bin_to_dec(binary_array)
[rows, cols] = size(binary_array);

% 检查行数是否是4的倍数
if mod(rows, 4) ~= 0
    error('行数必须是4的倍数');
end

% 计算结果数组的大小
new_rows = rows / 4;

% 创建一个 (n/4)*2 的十进制数组
decimal_array = zeros(new_rows, cols);

% 将二进制数组转换为十进制数组
for j = 1:cols
    for i = 1:new_rows
        % 将每4行的二进制数转换为字符串
        bin_str = num2str(binary_array(4*(i-1)+1:4*i, j)');
        % 去掉字符串中的空格
        bin_str = bin_str(bin_str ~= ' ');
        % 将二进制字符串转换为十进制数
        decimal_array(i, j) = bin2dec(bin_str);
    end
end
end

function [best_path,LLR] = kbest_soft(H, y, symbset, K,N0)

num = length(symbset);
nBitsPerSymbol = log2(num);
bit_mapping = qamdemod(symbset,16,'OutputType','bit','UnitAveragePower', true);
bit_mapping = bit_mapping';
[Q, R] = qr(H, 0);
y_tilde = Q' * y;
N = size(H,2);
S = symbset;
d = 0;

paths = zeros(K, N);    % 存储候选路径
path_metrics = inf(K, 1); % 路径度量 (初始化为无穷大)
path_metrics(1) = 0;    % 第一条路径度量初始化为 0
metrics_per_layer = zeros(K,N);

% 从底层到顶层递归
for i = N:-1:1
    candidates = []; % 临时存储扩展路径
    metrics = [];    % 对应的路径度量
    
    % 扩展每条路径
    for k = 1:K
        for s = S
            new_path = paths(k, :);
            new_path(i) = s; % 更新当前符号
            
            % 计算新的路径度量
            partial_metric = path_metrics(k) + abs(y_tilde(i) - R(i, i:end) * new_path(i:end).')^2;
            
            % 保存扩展的路径及度量
            candidates = [candidates; new_path];
            metrics = [metrics; partial_metric];
        end
    end
    
    % 对路径按度量排序，并保留前 K 条
    [sorted_metrics, idx] = sort(metrics);
    paths = candidates(idx(1:K), :);
    path_metrics = sorted_metrics(1:K);
    if i == N
        metrics_per_layer(:,i) = path_metrics;
    else 
        metrics_per_layer(:,i) = path_metrics -metrics_per_layer(:,i+1);
    end
end

% 输出最优路径及其度量
best_path = paths(1, :);
best_metric = path_metrics(1);

% 计算 LLR
% path_bit = qamdemod(paths,16,'OutputType','bit','UnitAveragePower', true);
% LLR = zeros(4, N);
LLR = inf(nBitsPerSymbol, N);

% 遍历每个符号位置（从第 1 层到第 N 层）
for i = 1:N
    % 获取当前层的候选路径和对应度量
    current_paths = paths(:, i); % 第 i 层的所有路径符号
    current_metrics = metrics_per_layer(:, i); % 对应的路径度量
    
    % 针对每个比特位置 (b = 1 到 nBitsPerSymbol)
    for b = 1:nBitsPerSymbol
        % 将符号路径映射到当前比特的值
        bits = zeros(size(current_paths));
        for k = 1:K
            symbol = current_paths(k);
            symbol_idx = find(S == symbol); % 找到符号在星座图中的索引
            bits(k) = bit_mapping(symbol_idx, b); % 映射到当前比特
        end
        
        % 计算当前比特为 0 和 1 的最小路径度量
        d0_candidates = current_metrics(bits == 0); % 比特为 0 的路径度量
        d1_candidates = current_metrics(bits == 1); % 比特为 1 的路径度量
        
        % 如果比特为 0 或比特为 1 的路径不存在，设置度量为无穷大
        if isempty(d0_candidates)
            d0 = 2/sqrt(10);
        else
            d0 = min(d0_candidates); % 最小路径度量
        end
        
        if isempty(d1_candidates)
            d1 = 2/sqrt(10);
        else
            d1 = min(d1_candidates); % 最小路径度量
        end
        % 计算 LLR
        LLR(b, i) = -1/N0*(d0 - d1);
    end
end

% 输出 LLR
% disp('Soft Output LLR:');
% disp(LLR);
% disp('Best Path:');
% disp(best_path);
% disp('Best Metric:');
% disp(best_metric);
    
end
