clear all;
close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
k = log2(sym_QAM);
iteration = 500;
num_symbol = 400;

trellis = poly2trellis(7,[171 133]);
tbl = 32;
rate = 1/2;

%EbN0s = 0:1:20;
SNR_dB = 0:5:20; % in dB
SNR = 10.^(SNR_dB./10);

qam_symbol =0:1:sym_QAM-1;
qam_signal = qammod(qam_symbol, sym_QAM, 'UnitAveragePower', true);
symbset_b3_0 = qam_signal(1:8);
symbset_b3_1 = qam_signal(9:16);
symbset_b2_0 = qam_signal([1:4,9:12]);
symbset_b2_1 = qam_signal([5:8,13:16]);
symbset_b1_0 = qam_signal([1,2,5,6,13,14,9,10]);
symbset_b1_1 = qam_signal([3,4,7,8,15,16,11,12]);
symbset_b0_0 = qam_signal([1,5,13,9,3,7,15,11]);
symbset_b0_1 = qam_signal([2,4,6,8,10,12,14,16]);

%y_QAM = zeros(1,iteration);
QAM_demod_SD = zeros(Ntx,num_symbol);
QAM_demod_SD_soft = zeros(Ntx,num_symbol);

errors_SD = zeros(length(SNR_dB),1);
errors_SD_soft = zeros(length(SNR_dB),1);
noise_variance = zeros(length(SNR_dB),1);

% X3 = zeros(Ntx,num_symbol*2);
% X2 = zeros(Ntx,num_symbol*2);
% X1 = zeros(Ntx,num_symbol*2);
% X0 = zeros(Ntx,num_symbol*2);


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

       
        % noise = sqrt(N0/2)*(randn(Nrx,num_symbol)+1i*randn(Nrx,num_symbol));
        % y = H * txsignal_QAM + noise;
        % noise = sqrt(N0/2)*(randn(Nrx,num_symbol)+1i*randn(Nrx,num_symbol));
        % y = H * txsignal_QAM' + noise; 

        y_withoutNoise = H * txsignal_QAM';
        y = awgn(y_withoutNoise,SNR_dB(l));

        w_MMSE = (H'*H+N0*eye(Ntx))\H';
        y_MMSE = w_MMSE * y;
        %SD
        near_Q = find_near_Q(y_MMSE,1/sqrt(10));
        [Q,R] = qr(H);
        radius_squared = norm(R *(near_Q-y_MMSE))^2;
        radius = sqrt(radius_squared);
        x_hat_SD_bit = zeros(Ntx, num_symbol*2);
        index_symbol = 1;
        while index_symbol <= num_symbol*2
            [r,count,lambda] = sphdec(H, y(:, index_symbol), symbset_b3_0, radius);
            lambda_c0 = lambda;
            [r,count,lambda] = sphdec(H, y(:, index_symbol), symbset_b3_1, radius);
            lambda_c1 = lambda;
            X3 = (lambda_c0-lambda_c1)/N0;

            [r,count,lambda] = sphdec(H, y(:, index_symbol), symbset_b2_0, radius);
            lambda_c0 = lambda;
            [r,count,lambda] = sphdec(H, y(:, index_symbol), symbset_b2_1, radius);
            lambda_c1 = lambda;
            X2 = (lambda_c0-lambda_c1)/N0;

            [r,count,lambda] = sphdec(H, y(:, index_symbol), symbset_b1_0, radius);
            lambda_c0 = lambda;
            [r,count,lambda] = sphdec(H, y(:, index_symbol), symbset_b1_1, radius);
            lambda_c1 = lambda;
            X1 = (lambda_c0-lambda_c1)/N0;

            [r,count,lambda] = sphdec(H, y(:, index_symbol), symbset_b0_0, radius);
            lambda_c0 = lambda;
            [r,count,lambda] = sphdec(H, y(:, index_symbol), symbset_b0_1, radius);
            lambda_c1 = lambda;
            X0 = (lambda_c0-lambda_c1)/N0;

            if r == 0
                radius = radius * sqrt(2);
            else
                x_hat_SD_bit(:, 4*index_symbol-3:4*index_symbol) = [X3,X2,X1,X0];
                index_symbol = index_symbol + 1;
            end
        end


        % soft = soft_decision_2x2(y_MMSE);
        % % SINR1 = (abs((w_MMSE(1,:)*H(:,1)))^2)/(abs((w_MMSE(1,:)*H(:,2)))^2+w_MMSE(1,:)*w_MMSE(1,:)'*N0);
        % % SINR2 = (abs((w_MMSE(2,:)*H(:,2)))^2)/(abs((w_MMSE(2,:)*H(:,1)))^2+w_MMSE(2,:)*w_MMSE(2,:)'*N0);
        % % QAM_demod_soft = soft .* [SINR1;SINR2];
        % QAM_demod_soft = soft;
        

        %解调
        % QAM_demod_SD_soft = qamdemod(x_hat_SD_bit',sym_QAM, 'OutputType','approxllr','UnitAveragePower', true, 'NoiseVariance', N0);
        % % QAM_demod_SD = soft_decode_2x2(y_MMSE');
        % QAM_demod_SD = qamdemod(x_hat_SD_bit',sym_QAM,'OutputType','bit','UnitAveragePower', true);

        QAM_demod_SD_soft = x_hat_SD_bit';
        QAM_demod_SD = zeros(size(QAM_demod_SD_soft));
        QAM_demod_SD(QAM_demod_SD_soft>0)=1;

        symbol_hard1 = vitdec(QAM_demod_SD(:,1),trellis,tbl,'cont','hard');
        symbol_hard2 = vitdec(QAM_demod_SD(:,2),trellis,tbl,'cont','hard');
        symbol_soft1 = vitdec(QAM_demod_SD_soft(:,1),trellis,tbl,'cont','unquant');
        symbol_soft2 = vitdec(QAM_demod_SD_soft(:,2),trellis,tbl,'cont','unquant');
        symbol_hard = [symbol_hard1,symbol_hard2];
        symbol_soft = [symbol_soft1,symbol_soft2];

        symbol_soft_dec = bin_to_dec(symbol_soft);
        symbol_hard_dec = bin_to_dec(symbol_hard);
        symbol_soft_dec = symbol_soft_dec(tbl/k+1:end);
        symbol_hard_dec = symbol_hard_dec(tbl/k+1:end);

        %计算误码
        errors_SD_soft(l) = errors_SD_soft(l) + sum(sum(symbol_soft_dec~=txsymbol_QAM(1:end-tbl/k)));
        errors_SD(l) = errors_SD(l) + sum(sum(symbol_hard_dec~=txsymbol_QAM(1:end-tbl/k)));
    end

end

%计算误码率
error_rate_soft = errors_SD/(iteration*num_symbol*Ntx);
error_rate_MMSE = errors_SD_soft/(iteration*num_symbol*Ntx);

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

function [demod_dec] = bin_to_dec2(x)
symbol_demod1 = x(1:4,:);
symbol_demod2 = x(5:8,:);
dec = [2^3;2^2;2^1;2^0];
symbol_demod_dec1 =  symbol_demod1.*dec ;
symbol_demod_dec2 =  symbol_demod2.*dec ;
demod_dec1 = sum(symbol_demod_dec1,1);
demod_dec2 = sum(symbol_demod_dec2,1);
demod_dec = [demod_dec1;demod_dec2];
end

function [bit] = llr_to_bit(x)
bit = x <0;
bit = double(bit);

end

function [x_hard] = demod_hard(x,d)
X_r=real(x); 
X_i=imag(x);
x_hard_r = zeros(size(x));
x_hard_i = zeros(size(x));

x_hard_r(X_r < -2*d) = -3;
x_hard_r(X_r > 2*d) = 3;
x_hard_r(X_r <= 0 & X_r >= -2*d) = -1;
x_hard_r(X_r > 0 & X_r <= 2*d) = 1;

x_hard_i(X_i < -2*d) = -3;
x_hard_i(X_i > 2*d) = 3;
x_hard_i(X_i <= 0 & X_i >= -2*d) = -1;
x_hard_i(X_i > 0 & X_i <= 2*d) = 1;

x_hard = x_hard_r + x_hard_i * 1i;
end

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

function [x4_soft] = soft_decode(x)
% x=txsignal_QAM(:).'; 
X_r=real(x); 
X_i=imag(x);

d = 1/sqrt(10);
% d = 1;
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
% X_int_bin = X >0;
% bin_to_dec = [2^3;2^2;2^1;2^0];
% X_int_dec = X_int_bin.' * bin_to_dec;
% x4_soft = reshape((X_int_dec(:)),2,[]); 
x4_soft = reshape((X(:)),2,[]);
x4_soft=  -x4_soft';
end


function [soft_2x2] = soft_decode_2x2(x)
[rows,cols] = size(x);
new_rows = 4 * rows;
soft_2x2 = zeros(new_rows,cols);
for i =1:rows
    soft_2x2(4*i-3:4*i,:) = soft_decode_1(x(i,:));
end

end

function [x4_soft] = soft_decode_1(x)

% x=txsignal_QAM(:).'; 
X_r=real(x); 
X_i=imag(x);

d = 1/sqrt(10);
% d = 1;
X_1 = zeros(size(x));
X_3 = zeros(size(x));% 初始化结果数组
for i = 1:length(x)
    if X_i(i) < -2*d
        X_1(i) = -2*d - 2*X_i(i);
    elseif X_i(i) > 2*d
        X_1(i) = 2*d - 2*X_i(i);
    else
        X_1(i) = -X_i(i);
    end
end

for i = 1:length(x)
    if X_r(i) < -2*d
        X_3(i) = 2*d + 2*X_r(i);
    elseif X_r(i) > 2*d
        X_3(i) = -2*d + 2*X_r(i);
    else
        X_3(i) = X_r(i);
    end
end

% if X_i < -2*d
%     X_1 = -2*d - 2*X_i;
% elseif X_i > 2*d
%     X_1 = 2*d - 2*X_i;
% else
%     X_1 = -X_i;
% end
% 
% if X_r < -2*d
%     X_3 = 2*d + 2*X_r;
% elseif X_r > 2*d
%     X_3 = -2*d + 2*X_r;
% else 
%     X_3 = X_r;
% end

X_2 = 2*d-abs(X_r);
X_0 = 2*d-abs(X_i);
X=[X_3; X_2;  X_1; X_0];  
X = X * 4*d;
% X_int_bin = X >0;
% bin_to_dec = [2^3;2^2;2^1;2^0];
% X_int_dec = X_int_bin.' * bin_to_dec;
% x4_soft = reshape((X_int_dec(:)),2,[]); 
% x4_soft = reshape((X(:)),2,[]);
x4_soft=  -X;
end


function [r,count, lambda] = sphdec(H, y, symbset, radius)

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
global RETVAL_SET;
SPHDEC_RADIUS = radius;
global lambda_layer2;
global lambda_layer1;

RETVAL        = zeros(n, 1);
TMPVAL        = zeros(n, 1);
SYMBSETSIZE   = length(symbset(:));
SEARCHFLAG    = 0;
count_SD  = 0;
sphdec_core(z, R, symbset, n, 0);



if SEARCHFLAG > 0
    r = RETVAL;
    lambda = [lambda_layer1;lambda_layer2];
    count = count_SD;
else
    r = 0;
    lambda = SPHDEC_RADIUS;
    count = count_SD;
end


clear SPHDEC_RADIUS RETVAL SYMBSETSIZE SEARCHFLAG count_SD lambda_layer2 lambda_layer1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sphdec_core(z, R, symbset, layer, dist)

global SPHDEC_RADIUS;
global RETVAL;
global TMPVAL;
global SYMBSETSIZE;
global SEARCHFLAG;
global count_SD;
global RETVAL_SET;
global lambda_layer2;
global lambda_layer1;

if (layer == 1)
    for ii = 1:SYMBSETSIZE
        TMPVAL(1) = symbset(ii);
        d = abs(z(1) - R(1,:)*TMPVAL)^2 + dist;
        count_SD = count_SD + 1;
        if (d <= SPHDEC_RADIUS)
            RETVAL        =  TMPVAL;
            SPHDEC_RADIUS =  d;
            SEARCHFLAG    =  SEARCHFLAG + 1;
            lambda_layer1 = d;
            RETVAL_SET = [RETVAL_SET,RETVAL];
        end
    end
else
    for ii = 1:SYMBSETSIZE
        TMPVAL(layer) = symbset(ii);
        d = abs(z(layer) - R(layer,[layer:end])*TMPVAL(layer:end))^2 + dist;
        count_SD = count_SD + 1;
        if (d <= SPHDEC_RADIUS)
            lambda_layer2 = d;
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
