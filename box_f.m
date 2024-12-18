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

txsymbol_QAM = [5;13]; % 生成0到15之间的整数作为符号
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
r = box_dec(H,txsignal_QAM,k_box);
demod = qamdemod(r,sym_QAM, 'UnitAveragePower', true);

function r = box_dec(H,y,k_box)

[Q, R] = qr(H, 0);
z = Q'*y;
n = size(H,2);
a = zeros(n,1);
s = zeros(n,k_box);

for layer = n:-1:1
    if layer ==n
        a(n) = z(layer)/R(layer,layer);
        s(n,:) = get_box(a(n));

    elseif layer ==1
        [row,col] = size(s);
        d_min = inf;
        for index = 1:col         
            
            a(1) = (z(1)-R(1,2)*s(layer+1,index))/R(1,1);
            s(1,index) = find_near_Q(a(1),1/sqrt(10));
            d = cal_distance(2,z,s(:,index),R);
            if d < d_min
                d_min = d;
                index_min = index;

            end
            r= s(:,index_min);

        end

        

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