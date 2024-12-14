clear;
close all;

M = 64; % 64QAM
box_num = 16;
n = sqrt(M); % 每一维的点数（假设M为正方形的QAM）
% 定义16-QAM的格雷码映射
% grayCode = ["0000", "0100", "1100", "1000";
%             "0001", "0101", "1101", "1001";
%             "0011", "0111", "1111", "1011";
%             "0010", "0110", "1110", "1010"];
% grayCode = ["0000", "1001", "0110", "1111";
%             "1010", "0011", "1100", "0101";
%             "0111", "1110", "0001", "1000";
%             "1101", "0100", "1011", "0010"];
% grayCode = ["0000", "0101", "1011", "1110";
%             "0011", "0110", "1100", "1001";
%             "1111", "1010", "0100", "0001";
%             "1101", "1000", "0010", "0111"];
new_map = find_map(M,box_num);
% new_map = [2    11     6    13    15     8     1     3     7     5    14     4    12     0    10     9];
grayCode = convertToBinaryMatrix(new_map);




% 定义QAM坐标
x = -(n-1):2:(n-1); % In-Phase (实部)
y = (n-1):-2:-(n-1); % Quadrature (虚部)

% 绘制16-QAM符号点
figure;
hold on;
for i = 1:n
    for j = 1:n
        % 绘制星座点
        plot(x(j), y(i), 'b*', 'MarkerSize', 10);
        % 在星座点旁标注格雷码
        text(x(j)+0.2, y(i), grayCode(i, j), 'HorizontalAlignment', 'left', 'FontSize', 10);
        % 在星座点下方标注索引值
        index = binaryToDecimal(grayCode(i,j)); % 索引从0开始
        text(x(j)-0.2, y(i)-0.3, num2str(index), 'HorizontalAlignment', 'right', 'FontSize', 10);
    end
end

% 设置图形标题和轴标签
title([num2str(M) '-QAM Symbol Mapping']);
xlabel('In-Phase');
ylabel('Quadrature');
axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1]); % 动态设置坐标范围
grid on;
hold off;


function decimal = binaryToDecimal(binaryStr)
    % 将二进制字符串转换为十进制数
    decimal = bin2dec(binaryStr);
end

function binaryMatrix = convertToBinaryMatrix(new_map)
    n = sqrt(length(new_map));
    
    % 初始化 4x4 矩阵
    binaryMatrix = strings(n, n);
    
    % 转换为二进制字符串并填充矩阵
    for i = 1:n^2
        binaryMatrix(i) = dec2bin(new_map(i), log2(n^2)); % 转换为4位二进制
    end
    
    % 重塑为 4x4 矩阵
    binaryMatrix = reshape(binaryMatrix, [n, n]);
end

function return_bit = find_map(M,box_num)
% 参数设置
num_trials = 1e6; % 随机尝试次数
n = sqrt(box_num); % 正方形边长，改为3x3的矩形区域（面积为9）
bit_num = log2(M);

% 定义矩形的索引
rows = sqrt(M);
cols = rows;
rectangles = [];
for i = 1:rows-n+1
    for j = 1:cols-n+1
        rect = [];
        for di = 0:n-1
            for dj = 0:n-1
                rect = [rect, (i+di-1)*cols + (j+dj)];
            end
        end
        rectangles = [rectangles; rect];
    end
end

success = false; % 标志位
for trial = 1:num_trials
    % 随机生成比特映射
    bits = randperm(M) - 1; % 生成0到15的随机排列
    % bits = [9    10     5     1     8    15     3     7    12    14    13     6     2     0    11     4];
    
    % 检查每个矩形的条件
    is_valid = true;
    for r = 1:size(rectangles, 1)
        rect = rectangles(r, :); % 当前矩形的点
        bits_in_rect = zeros(1,bit_num);
        for bit = 1:bit_num
            for rect_num = 1:length(rect)
                rrt = bits(rect(rect_num));
                bits_in_rect(rect_num) = bitget(bits(rect(rect_num)), bit);
            end
            if all(bits_in_rect == 0) || all(bits_in_rect == 1)
                is_valid = false;
                break;
            end
        end
    end
    
    if is_valid
        success = true;
        return_bit = bits;
        break;
    end
end
end
