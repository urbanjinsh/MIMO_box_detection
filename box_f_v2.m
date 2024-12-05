clear;
close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
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

txsymbol_QAM = [1;9]; % 生成0到15之间的整数作为符号
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
[r,count] = box_dec(H,txsignal_QAM,k_box);
demod = qamdemod(r,sym_QAM, 'UnitAveragePower', true);

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

