clear; close all

sym_QAM = 16;  
k = log2(sym_QAM);
iteration = 400;
num_symbol = 200;

trellis = poly2trellis(7,[171 133]);
tbl = 32;
rate = 1/2;

Ntx = 2;
Nrx = 2;


array = randi([0,sym_QAM-1], Ntx, num_symbol);
% 获取数组的大小
[rows, cols] = size(array);

% 将数组中的每个元素转换为4位二进制字符串
binary_strings = dec2bin(array, 4);

% 创建一个2*4n的double型数组
binary_array = zeros(rows, 4 * cols);

% 将二进制字符串转换为数字，并存储在2*4n的double型数组中
for i = 1:rows
    for j = 1:cols
        % 获取当前元素的二进制字符串
        bin_str = binary_strings((i-1)*cols + j, :);
        % 将二进制字符串的每一位转换为数字，并存储在结果数组中
        for k = 1:4
            binary_array(i, 4*(j-1) + k) = str2double(bin_str(k));
        end
    end
end

txsymbol_QAM_enc = convenc(binary_array,trellis);
binary_array2=[1,0,1,1];
% txsymbol_QAM_enc1 = convenc(binary_array2,trellis);
