clear all;
close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
iteration = 400;
num_symbol = 200;
k_best_num = 4;
k_box_num = 4;

SNR_dB = 0:4:20; % in dB
SNR = 10.^(SNR_dB./10);

errors_KB = zeros(length(SNR_dB),1);
errors_MMSE = zeros(length(SNR_dB),1);
errors_box = zeros(length(SNR_dB),1);

qam_symbol =0:1:sym_QAM-1;
qam_signal = qammod(qam_symbol, sym_QAM, 'UnitAveragePower', true); 

PED_count_box = zeros(length(SNR_dB),1);
PED_count_KB = zeros(length(SNR_dB),1);


for l = 1:length(SNR_dB)
    N0 = 1/(10^(SNR_dB(l)/10));
    for iter = 1:iteration
        txsymbol_QAM = randi([0,sym_QAM-1], Ntx, num_symbol); % 生成0到15之间的整数作为符号
        txsignal_QAM = qammod(txsymbol_QAM, sym_QAM, 'UnitAveragePower', true); % 使用灰度映射调制符号
        H = sqrt(0.5)*(randn(Nrx,Ntx)+1i*randn(Nrx,Ntx));
        noise = sqrt(N0/2)*(randn(Nrx,num_symbol)+1i*randn(Nrx,num_symbol));
        % noise = 0;
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
        
        %box
        x_hat_box = zeros(Ntx, num_symbol);
        index_symbol = 1;
        while index_symbol <= num_symbol
            [r,count] = box_dec(H, y(:, index_symbol), k_box_num);
            
            x_hat_box(:, index_symbol) = r;
            PED_count_box(l) = PED_count_box(l) + count;
            index_symbol = index_symbol + 1;
        end


        QAM_demod_MMSE = qamdemod(y_MMSE,sym_QAM, 'UnitAveragePower', true);
        QAM_demod_KB = qamdemod(x_hat_KB,sym_QAM, 'UnitAveragePower', true);
        QAM_demod_box = qamdemod(x_hat_box,sym_QAM, 'UnitAveragePower', true);
        errors_MMSE(l) = errors_MMSE(l) + sum(sum(QAM_demod_MMSE~=txsymbol_QAM));
        errors_KB(l) = errors_KB(l) + sum(sum(QAM_demod_KB~=txsymbol_QAM));
        errors_box(l) = errors_box(l) + sum(sum(QAM_demod_box~=txsymbol_QAM));
        
    end
    disp(l)
end

error_rate_MMSE = errors_MMSE/(iteration*num_symbol*Ntx);
error_rate_box = errors_box/(iteration*num_symbol*Ntx);
error_rate_KB = errors_KB/(iteration*num_symbol*Ntx);

figure;
semilogy(SNR_dB,error_rate_MMSE,'b-s',SNR_dB,error_rate_box,'g-x',SNR_dB,error_rate_KB,'r-o');
legend('MMSE','box,k=4','KB,k=4');
title('BER of  MMSE, box and KB');
xlabel('SNR in dB');
ylabel('BER');

figure; % 创建新图形窗口
plot(SNR_dB, PED_count_box, '-o'); % 绘制第一组数据
hold on; % 保持当前图形
plot(SNR_dB, PED_count_KB, '-x'); % 绘制第二组数据

% 添加图例
legend('box,k=4','KB,k=4');

% 添加标题和轴标签
title('PED times on different SNR');
xlabel('SNR (dB)');
ylabel('PED calculation tiumes');

% 显示网格
grid on;


function [r,count_KB] = kbest(H, y, symbset, k)

[Q, R] = qr(H, 0);
z = Q'*y;
n = size(H,2);

d = 0;
SYMBSETSIZE = length(symbset(:));
distance = zeros(n,k);
temp_distance = zeros(k,SYMBSETSIZE);
RETVAL        = zeros(n, 1);
TMPVAL        = zeros(n, 1);
TMPSET        = zeros(n, k);
k_best_symbset = zeros(n, k);
count_KB = 0;

for layer = n : -1 :1
    if layer == n
        for ii = 1:SYMBSETSIZE
            TMPVAL(layer) = symbset(ii);
            temp_distance(k,ii) = abs(z(layer) - R(layer,[layer:end])*TMPVAL(layer:end))^2 + temp_distance(k,ii);
            count_KB = count_KB+ 1;
        end
    
        [temp_distance_asc,sort_index] = sort(temp_distance(k,:));
        k_best_symbset(n,:) = symbset(sort_index(1:k));
        distance(n,:) = temp_distance(k,sort_index(1:k));

    else
        for i_kbest = 1:k
            TMPVAL(layer+1:end) = k_best_symbset(layer+1:end,i_kbest);
            for ii = 1:SYMBSETSIZE

                TMPVAL(layer) = symbset(ii);
                temp_distance(i_kbest,ii) = abs(z(layer) - R(layer,[layer:end])*TMPVAL(layer:end))^2 + distance(layer+1,i_kbest);
                count_KB = count_KB+ 1;
            end
            
        end
    
        [temp_distance_asc,sort_index] = sort(temp_distance(:));
        pos2 = mod(sort_index(1:k),4);
        pos2(pos2 == 0) = pos2(pos2 == 0) + 4;
        k_best_symbset(layer+1,:) = k_best_symbset(layer+1,pos2);
        pos2_layer1 = fix((sort_index(1:k)-0.1)/4)+1;
        k_best_symbset(layer,:) = symbset(pos2_layer1);
        
        r = k_best_symbset(:,1);
        
    end


    
end


end

function [best_path,count] = box_dec(H,y,k_box)
count = 0;
[Q, R] = qr(H, 0);
z = Q'*y;
n = size(H,2);
a = zeros(n,1);
s = zeros(n,k_box);
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
            partial_candidates = find_near_Q(a(1),1/sqrt(10));
            for s = partial_candidates
                new_path = paths(k, :);
                new_path(layer) = s; % 更新当前符号
                partial_metric = path_metrics(k) + abs(z(layer) - R(layer, layer:end) * new_path(layer:end).')^2;
                count = count+1;
                metrics = [metrics; partial_metric];
                candidates = [candidates; new_path];
            end      
        end
        metrics_per_layer(:,layer) = path_metrics -metrics_per_layer(:,layer+1);
        [sorted_metrics, idx] = sort(metrics);
        paths = candidates(idx(1:K), :);
        path_metrics = sorted_metrics(1:K);

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
best_path = paths(1, :);
best_metric = path_metrics(1);

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
    if abs(x_re) >= 3
        x_re = sign(x_re) * 2.5;
    end
    re1 = (fix(x_re/(2))-0.5)*2;
    re2 = (fix(x_re/(2))-0.5 + 1)*2;

    % 处理虚部
    if abs(x_im) >= 3
        x_im = sign(x_im) * 2.5;
    end
    im1 = (fix(x_im/(2))-0.5)*2;
    im2 = (fix(x_im/(2))-0.5 + 1)*2;

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

function [d,count] = cal_distance(layer,x,s,R)
d = 0; % 初始化距离
count = 0;
for i = layer:-1:1 % 从 N 层递归到第 1 层
    inner_sum = 0; % 初始化内层求和
    for j = i:layer
        inner_sum = inner_sum + R(i, j) * s(j); % 内层求和
        count = count +1;
    end
    d = d + abs(x(i) - inner_sum)^2; % 计算当前层的距离并累计
end
end