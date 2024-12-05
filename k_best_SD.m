clear all;
close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
iteration = 400;
num_symbol = 200;
k_best_num = 4;

SNR_dB = 0:4:20; % in dB
SNR = 10.^(SNR_dB./10);

errors_KB = zeros(length(SNR_dB),1);
errors_MMSE = zeros(length(SNR_dB),1);
errors_SD = zeros(length(SNR_dB),1);

qam_symbol =0:1:sym_QAM-1;
qam_signal = qammod(qam_symbol, sym_QAM, 'UnitAveragePower', true); 

PED_count_SD = zeros(length(SNR_dB),1);
PED_count_KB = zeros(length(SNR_dB),1);


for l = 1:length(SNR_dB)
    N0 = 1/(10^(SNR_dB(l)/10));
    for iter = 1:iteration
        txsymbol_QAM = randi([0,sym_QAM-1], Ntx, num_symbol); % 生成0到15之间的整数作为符号
        txsignal_QAM = qammod(txsymbol_QAM, sym_QAM, 'UnitAveragePower', true); % 使用灰度映射调制符号
        H = sqrt(0.5)*(randn(Nrx,Ntx)+1i*randn(Nrx,Ntx));
        noise = sqrt(N0/2)*(randn(Nrx,num_symbol)+1i*randn(Nrx,num_symbol));
        y = H * txsignal_QAM + noise;

        %MMSE
        w_MMSE = (H'*H+N0*eye(Ntx))\H';
        y_MMSE = w_MMSE * y;

        %K-Best
        x_hat_KB = zeros(Ntx, num_symbol);
        index_symbol = 1;
        while index_symbol <= num_symbol
            [r,count] = kbest(H, y(:, index_symbol), qam_signal, k_best_num);
            
            x_hat_KB(:, index_symbol) = r;
            PED_count_KB(l) = PED_count_KB(l) + count;
            index_symbol = index_symbol + 1;
        end
        
        %SD
        near_Q = find_near_Q(y_MMSE,1/sqrt(10));
        % Q is the Constellation diagram
        [Q,R] = qr(H);
        radius_squared = norm(R *(near_Q-y_MMSE))^2;
        radius = sqrt(radius_squared);
        x_hat_SD = zeros(Ntx, num_symbol);
        index_symbol = 1;
        while index_symbol <= num_symbol
            [r,count] = sphdec(H, y(:, index_symbol), qam_signal, radius);
            if r == 0
                radius = radius * sqrt(2);
            else
                x_hat_SD(:, index_symbol) = r;
                PED_count_SD(l) = PED_count_SD(l) + count;
                index_symbol = index_symbol + 1;
            end
        end

        QAM_demod_MMSE = qamdemod(y_MMSE,sym_QAM, 'UnitAveragePower', true);
        QAM_demod_KB = qamdemod(x_hat_KB,sym_QAM, 'UnitAveragePower', true);
        QAM_demod_SD = qamdemod(x_hat_SD,sym_QAM, 'UnitAveragePower', true);
        errors_MMSE(l) = errors_MMSE(l) + sum(sum(QAM_demod_MMSE~=txsymbol_QAM));
        errors_KB(l) = errors_KB(l) + sum(sum(QAM_demod_KB~=txsymbol_QAM));
        errors_SD(l) = errors_SD(l) + sum(sum(QAM_demod_SD~=txsymbol_QAM));
        
    end
    disp(l)
end

error_rate_MMSE = errors_MMSE/(iteration*num_symbol*Ntx);
error_rate_SD = errors_SD/(iteration*num_symbol*Ntx);
error_rate_KB = errors_KB/(iteration*num_symbol*Ntx);

figure;
semilogy(SNR_dB,error_rate_MMSE,'b-s',SNR_dB,error_rate_SD,'g-x',SNR_dB,error_rate_KB,'r-o');
legend('MMSE','SD','KB');
title('BER of  MMSE, SD and KB');
xlabel('SNR in dB');
ylabel('BER');

figure; % 创建新图形窗口
plot(SNR_dB, PED_count_SD, '-o'); % 绘制第一组数据
hold on; % 保持当前图形
plot(SNR_dB, PED_count_KB, '-x'); % 绘制第二组数据

% 添加图例
legend('SD','KB');

% 添加标题和轴标签
title('PED times on different SNR');
xlabel('SNR (dB)');
ylabel('PED calculation tiumes');

% 显示网格
grid on;

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

function [x_hard] = find_near_Q(x,d)
X_r=real(x); 
X_i=imag(x);
x_hard_r = zeros(size(x));
x_hard_i = zeros(size(x));

x_hard_r(X_r < -2*d) = -3*d;
x_hard_r(X_r > 2*d) = 3*d;
x_hard_r(X_r <= 0 & X_r >= -2*d) = -1*d;
x_hard_r(X_r > 0 & X_r <= 2*d) = 1*d;

x_hard_i(X_i < -2*d) = -3*d;
x_hard_i(X_i > 2*d) = 3*d;
x_hard_i(X_i <= 0 & X_i >= -2*d) = -1*d;
x_hard_i(X_i > 0 & X_i <= 2*d) = 1*d;

x_hard = x_hard_r + x_hard_i * 1i;

end

function [best_path,count] = kbest(H, y, symbset, K)
count = 0;
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
    if i ==N
        for s = S
            new_path = zeros(1, N);
            new_path(i) = s; % 更新当前符号
            
            % 计算新的路径度量
            partial_metric = abs(y_tilde(i) - R(i, i:end) * new_path(i:end).')^2;
            count = count+1;
            % 保存扩展的路径及度量
            candidates = [candidates; new_path];
            metrics = [metrics; partial_metric];
        end
        [sorted_metrics, idx] = sort(metrics);
        paths = candidates(idx(1:K), :); % 初始化前 K 条路径
        path_metrics = sorted_metrics(1:K);
    else
        for k = 1:K
            for s = S
                new_path = paths(k, :);
                new_path(i) = s; % 更新当前符号
            
                % 计算新的路径度量
                partial_metric = path_metrics(k) + abs(y_tilde(i) - R(i, i:end) * new_path(i:end).')^2;
                count = count+1;
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
end

% 输出最优路径及其度量
best_path = paths(1, :);
best_metric = path_metrics(1);

end


    
