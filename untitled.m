A = [6,5,4,3]; % 二进制表示为 1111
% bit = bitget(A, 1:6); % 获取最低位，结果为 1
bit1 = bitget(A, 1); % 获取最低位，结果为 1
bit2 = bitget(A, 2); % 获取第二位，结果为 1
bit3 = bitget(A, 3); % 获取第三位，结果为 1
bit4 = bitget(A, 4); % 获取第四位，结果为 1
bit5 = bitget(A, 5); % 获取第四位，结果为 1
bit6 = bitget(A, 6); % 获取第四位，结果为 1
bit = zeros(length(A), 6); % 初始化一个矩阵来存储结果

for i = 1:length(A)
    for j = 1:6
        bit(i, j) = bitget(A(i), j);
    end
end