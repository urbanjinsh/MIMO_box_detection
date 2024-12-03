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
r = box_dec(H,txsignal_QAM,k_box,N0);
demod = qamdemod(r,sym_QAM, 'UnitAveragePower', true);

function r = box_dec(H,y,k_box,n0)
sym_QAM = 16;
bit_num = log2(sym_QAM);
[Q, R] = qr(H, 0);
z = Q'*y;
z = [0.237639912667480 - 0.677230835384780i;-0.588410654381244 - 0.451478493395317i];
n = size(H,2);
a = zeros(n,1);
s = zeros(n,k_box);
s_bit = zeros(n,k_box*bit_num);
s_llr = zeros(n,bit_num);
l = zeros(n,k_box);

for layer = n:-1:1
    if layer ==n
        
        a(n) = z(layer)/R(layer,layer);
        
        s(n,:) = get_box(a(n));
        temp = qamdemod(s(n,:),sym_QAM,'OutputType','bit', 'UnitAveragePower', true);
        s_bit(n,:) = temp(:)';
        

        %for 循环每个bit
        for i = 1:bit_num
            lambda_c0_min = inf;
            lambda_c1_min = inf;
            bit_temp = s_bit(n,[i,i+bit_num,i+bit_num*2,i+bit_num*3]);
            [bit_temp_asc,sort_index] = sort(bit_temp);
            if bit_temp_asc(1) == bit_temp_asc(bit_num)
                if bit_temp_asc(1) == 1
                    s_llr(n,i) =  1/n0 *(10-0);
                
                elseif bit_temp_asc(1) == 0
                    s_llr(n,i) =  1/n0 *(0-10);
                end

            else
                for box_point = 1:4
                    if bit_temp(box_point) == 1
                        temp = [0;s(n,box_point)];
                        lambda_c1 = abs(z(n) - R(n,[n:end])*temp(n:end))^2;
                        if lambda_c1 < lambda_c1_min
                            lambda_c1_min = lambda_c1;
                            lambda_c1_min_idx = box_point;
                        end
                
                    elseif bit_temp(box_point) == 0
                        temp = [0;s(n,box_point)];
                        lambda_c0 = abs(z(n) - R(n,[n:end])*temp(n:end))^2;
                        if lambda_c0 < lambda_c0_min
                            lambda_c0_min = lambda_c0;
                            lambda_c0_min_idx = box_point;
                        end
                        
                    end
                end
                s_llr(n,i) =  1/n0 *(lambda_c0_min - lambda_c1_min);



            end

        end



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