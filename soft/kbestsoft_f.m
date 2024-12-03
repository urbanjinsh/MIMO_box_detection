clear;
close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;
bit_num = log2(16);
iteration = 400;
num_symbol = 1;

SNR_dB = 0:4:20; % in dB
SNR = 10.^(SNR_dB./10);

errors_ML = zeros(length(SNR_dB),1);
errors_MMSE = zeros(length(SNR_dB),1);
errors_SD = zeros(length(SNR_dB),1);

qam_symbol =0:1:sym_QAM-1;
qam_signal = qammod(qam_symbol, sym_QAM, 'UnitAveragePower', true); 

PED_count_SD = zeros(length(SNR_dB),1);
PED_count_ML = zeros(length(SNR_dB),1);

k_box = 4;
k_best_num = 4;

txsymbol_QAM = [6;1]; % 生成0到15之间的整数作为符号
txsignal_QAM = qammod(txsymbol_QAM, sym_QAM, 'UnitAveragePower', true); % 使用灰度映射调制符号
% H = [0.7,0.2;0.8,0.3];
H=[1.4 - 0.6i,0.7 - 0.7i;-0.8 - 0.6i,0.3 + 0.06i];
% H = [1,0;0,1];
l=6;
N0 = 1/(10^(SNR_dB(l)/10));
noise = sqrt(N0/2)*(randn(Nrx,num_symbol)+1i*randn(Nrx,num_symbol));
% noise = 0;
txsignal_QAM = H*txsignal_QAM+noise;


% box_list= get_box(txsignal_QAM(2));
r = kbest_soft(H,txsignal_QAM,qam_signal,k_best_num,N0);
demod = qamdemod(r,sym_QAM, 'UnitAveragePower', true);

function LLR = kbest_soft(H, y, symbset, K,N0)

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
            d0 = 100*N0;
        else
            d0 = min(d0_candidates); % 最小路径度量
        end
        
        if isempty(d1_candidates)
            d1 = 100*N0;
        else
            d1 = min(d1_candidates); % 最小路径度量
        end
        % 计算 LLR
        LLR(b, i) = 1/N0*(d0 - d1);
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