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

txsymbol_QAM = [3;13]; % 生成0到15之间的整数作为符号
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
[r,LLR] = box_dec(H,txsignal_QAM,k_box,N0,qam_signal);
demod = qamdemod(r,sym_QAM, 'UnitAveragePower', true);

function [best_path,LLR] = box_dec(H,y,k_box,N0,symbset)
count = 0;
[Q, R] = qr(H, 0);
z = Q'*y;
n = size(H,2);
a = zeros(n,1);
K = k_box;
N=n;
paths = zeros(K, N);    % 存储候选路径
path_metrics = inf(K, 1); % 路径度量 (初始化为无穷大)
path_metrics(1) = 0;    % 第一条路径度量初始化为 0
metrics_per_layer = zeros(K,N);

for layer = n:-1:1
    candidates = []; % 临时存储扩展路径
    metrics = [];    % 对应的路径度量
    if layer ==n
        a(n) = z(layer)/R(layer,layer);
        partial_candidates = get_box(a(n));
        for s = partial_candidates
            new_path = zeros(1, N);
            new_path(layer) = s; % 更新当前符号
            partial_metric = abs(z(layer) - R(layer, layer:end) * new_path(layer:end).')^2;
            count = count+1;
            metrics = [metrics; partial_metric];
            candidates = [candidates; new_path];
        end
        paths = candidates;
        path_metrics = metrics;
        metrics_per_layer(:,layer) = path_metrics;

    elseif layer ==1
        
        for k = 1:K
            a(1) = (z(layer) - R(layer, layer+1:N) * paths(k,layer+1:N)) / R(layer, layer);
            partial_candidates =  find_near_Q(a(1),1/sqrt(10));
            for s = partial_candidates
                new_path = paths(k, :);
                new_path(layer) = s; % 更新当前符号
                partial_metric = path_metrics(k) + abs(z(layer) - R(layer, layer:end) * new_path(layer:end).')^2;
                count = count+1;
                metrics = [metrics; partial_metric];
                candidates = [candidates; new_path];
            end      
        end
        
       
        paths = candidates;
        path_metrics = metrics;
        metrics_per_layer(:,layer) = path_metrics -metrics_per_layer(:,layer+1);

    else% layer !=n and 1
        
        for k = 1:K
            a(layer) = (x(layer) - r(layer, layer+1:N) * paths(k,layer+1:N)) / r(layer, layer);
            partial_candidates = get_box(a(n));
            for s = partial_candidates
                new_path = paths(k, :);
                new_path(layer) = s; % 更新当前符号
                partial_metric = path_metrics(k) + abs(z(layer) - R(layer, layer:end) * new_path(layer:end).')^2;
                count = count+1;
                metrics = [metrics; partial_metric];
                candidates = [candidates; new_path];
            end      
        end
        paths = candidates;
        path_metrics = metrics;
        metrics_per_layer(:,layer) = path_metrics -metrics_per_layer(:,layer+1);
    end


end
[best_metric,min_index] = min(path_metrics);
best_path = paths(min_index, :);


num = 16;
nBitsPerSymbol = log2(num);
S = symbset;
bit_mapping = qamdemod(symbset,16,'OutputType','bit','UnitAveragePower', true);
bit_mapping = bit_mapping';
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
            d0 = 1/sqrt(10);
        else
            d0 = min(d0_candidates); % 最小路径度量
        end
        
        if isempty(d1_candidates)
            d1 = 1/sqrt(10);
        else
            d1 = min(d1_candidates); % 最小路径度量
        end
        % 计算 LLR
        LLR(b, i) = -1/N0*(d0 - d1);
    end
end

end

function box_list = get_box(x) %% equation 10 ??
d = 1/sqrt(10);
n=1;
x=x/d;
x_re = real(x);
x_im = imag(x);

box_list = zeros(4 * n, 1);

for i = 1:n
    % 获取当前复数的实部和虚部
    x_re = real(x(i));
    x_im = imag(x(i));

    % 处理实部
    re_list = get_near_two(x_re);
    re1 = re_list(1);
    re2 = re_list(2);

    % 处理虚部
    im_list = get_near_two(x_im);
    im1 = im_list(1);
    im2 = im_list(2);

    % 生成4个不同的复平面点
    box_list(4 * (i - 1) + 1) = re1 + 1i * im1;
    box_list(4 * (i - 1) + 2) = re1 + 1i * im2;
    box_list(4 * (i - 1) + 3) = re2 + 1i * im1;
    box_list(4 * (i - 1) + 4) = re2 + 1i * im2;
    
end
box_list = box_list.';
box_list = box_list * d;
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

function d = cal_distance(layer,x,s,R)
d = 0; % 初始化距离
% for i = layer:-1:1 % 从 N 层递归到第 1 层
%     inner_sum = 0; % 初始化内层求和
%     for j = i:layer
%         inner_sum = inner_sum + R(i, j) * s(j); % 内层求和
%     end
%     d = d + abs(x(i) - inner_sum)^2; % 计算当前层的距离并累计
% end
for i = layer:-1:1 % 从 N 层递归到第 1 层
    d= abs(x(i) - R(i,[i:end])*s(i:end))^2 + d;
end

end

function near_list = get_near_two(x)
near_list=zeros(1,2);
if x <=-1
    near_list = [-3,-1];
elseif x<1 && x > -1
    near_list = [-1,1];
elseif x >=1
    near_list = [1,3];

end
end