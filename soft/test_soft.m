clear; close all

sym_QAM = 16;  
k = log2(sym_QAM);
iteration = 400;
num_symbol = 200;
Ntx = 2;
txsymbol_QAM = randi([0,sym_QAM-1],num_symbol,Ntx);
H=[1,0.5;0.3,1];

symbol = [3,2,5;14,2,5]; 
% symbol = [0;0;1;0];
signal = qammod(txsymbol_QAM,sym_QAM,'gray');
signal = awgn(signal,6);
% channel_signal = signal*H  ;
% 
% w_ZF = (H'*H)\H';
% orignal_signal =  channel_signal *w_ZF;

signal_soft = qamdemod(signal,sym_QAM,'OutputType','approxllr','NoiseVariance', 1);
signal_soft_own = soft_decode_2x2(signal);




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

% d = 1/sqrt(10);
d = 1;
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