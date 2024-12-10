clear;
close all;

% 定义16-QAM的格雷码映射
grayCode = ["0000", "0100", "1100", "1000";
            "0001", "0101", "1101", "1001";
            "0011", "0111", "1111", "1011";
            "0010", "0110", "1110", "1010"];
grayCode = ["0000", "1001", "0110", "1111";
            "1010", "0011", "1100", "0101";
            "0111", "1110", "0001", "1000";
            "1101", "0100", "1011", "0010"];
grayCode = ["0000", "0101", "1011", "1110";
            "0011", "0110", "1100", "1001";
            "1111", "1010", "0100", "0001";
            "1101", "1000", "0010", "0111"];
new_map = [2,15,4,8,10,6,3,14,1,7, 0,9,11 , 12 ,   13 ,   5];
grayCode = convertToBinaryMatrix(new_map);



% 定义QAM坐标
x = [-3, -1, 1, 3]; % In-Phase (实部)
y = [3, 1, -1, -3]; % Quadrature (虚部)

% 绘制16-QAM符号点
figure;
hold on;
for i = 1:4
    for j = 1:4
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
title('16-QAM Symbol Mapping');
xlabel('In-Phase');
ylabel('Quadrature');
axis([-4 4 -4 4]); % 设置坐标范围
grid on;
hold off;


function decimal = binaryToDecimal(binaryStr)
    % 将二进制字符串转换为十进制数
    decimal = bin2dec(binaryStr);
end

function binaryMatrix = convertToBinaryMatrix(new_map)
    % 检查输入是否为长度为16的数组
    if length(new_map) ~= 16
        error('Input array must have exactly 16 elements.');
    end
    
    % 初始化 4x4 矩阵
    binaryMatrix = strings(4, 4);
    
    % 转换为二进制字符串并填充矩阵
    for i = 1:16
        binaryMatrix(i) = dec2bin(new_map(i), 4); % 转换为4位二进制
    end
    
    % 重塑为 4x4 矩阵
    binaryMatrix = reshape(binaryMatrix, [4, 4]);
end