clear;
close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
iteration = 400;
num_symbol = 200;

SNR_dB = 0:4:20; % in dB
SNR = 10.^(SNR_dB./10);

errors_ML = zeros(length(SNR_dB),1);
errors_MMSE = zeros(length(SNR_dB),1);
errors_SD = zeros(length(SNR_dB),1);

qam_symbol =0:1:sym_QAM-1;
qam_signal = qammod(qam_symbol, sym_QAM, 'UnitAveragePower', true); 

PED_count_SD = zeros(length(SNR_dB),1);
PED_count_ML = zeros(length(SNR_dB),1);

k_best_num = 4;

txsymbol_QAM = [8;6]; % 生成0到15之间的整数作为符号
txsignal_QAM = qammod(txsymbol_QAM, sym_QAM, 'UnitAveragePower', true); % 使用灰度映射调制符号
H = [0.7,0.3;0.3,0.78];
% H = [1,0;0,1];

function [r,count] = sphdec(H, y, symbset, radius)

if nargin == 3
    radius = realmax;
end

if size(H, 1) < size(H, 2)
	H = [H; zeros(size(H, 2) - size(H, 1), size(H, 2))];
end

[Q, R] = qr(H, 0);
z = Q'*y;
n = size(H,2);


% add examine this variable before make it global
global SPHDEC_RADIUS;
global RETVAL;
global TMPVAL;
global SYMBSETSIZE;
global SEARCHFLAG;
global count_SD;
SPHDEC_RADIUS = radius;

RETVAL        = zeros(n, 1);
TMPVAL        = zeros(n, 1);
SYMBSETSIZE   = length(symbset(:));
SEARCHFLAG    = 0;
count_SD  = 0;
sphdec_core(z, R, symbset, n, 0);

if SEARCHFLAG > 0
    r = RETVAL;
    count = count_SD;
else
    r = 0;
    count = count_SD;
end


clear SPHDEC_RADIUS RETVAL SYMBSETSIZE SEARCHFLAG count_SD;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sphdec_core(z, R, symbset, layer, dist)

global SPHDEC_RADIUS;
global RETVAL;
global TMPVAL;
global SYMBSETSIZE;
global SEARCHFLAG;
global count_SD;

if (layer == 1)
    for ii = 1:SYMBSETSIZE
        TMPVAL(1) = symbset(ii);
        d = abs(z(1) - R(1,:)*TMPVAL)^2 + dist;
        count_SD = count_SD + 1;
        if (d <= SPHDEC_RADIUS)
            RETVAL        =  TMPVAL;
            SPHDEC_RADIUS =  d;
            SEARCHFLAG    =  SEARCHFLAG + 1;
        end
    end
else
    for ii = 1:SYMBSETSIZE
        TMPVAL(layer) = symbset(ii);
        d = abs(z(layer) - R(layer,[layer:end])*TMPVAL(layer:end))^2 + dist;
        count_SD = count_SD + 1;
        if (d <= SPHDEC_RADIUS)
            sphdec_core(z, R, symbset, layer-1, d);
        end
    end
end




end



r = kbest(H,txsignal_QAM,qam_signal,k_best_num);

demod = qamdemod(r, sym_QAM, 'UnitAveragePower', true);


function best_path = kbest(H, y, symbset, K)

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

end




