clear;
close all;


% 参数设置
M = 64; % 64QAM
num_trials = 1e6; % 随机尝试次数
bit_num = log2(M);

% 定义矩形的索引
rows = sqrt(M);
cols = rows;
rectangles = [];
for i = 1:rows-1
    for j = 1:cols-1
        idx = [(i-1)*cols + j, (i-1)*cols + j+1, ...
               i*cols + j, i*cols + j+1];
        rectangles = [rectangles; idx];
    end
end

success = false; % 标志位
for trial = 1:num_trials
    % 随机生成比特映射
    bits = randperm(M) - 1; % 生成0到15的随机排列
    
    % 检查每个矩形的条件
    is_valid = true;
    for r = 1:size(rectangles, 1)
        rect = rectangles(r, :); % 当前矩形的点
        bits_in_rect = zeros(1,bit_num);
        for bit = 1:length(rect)
            for rect_num = 1:bit_num
                bits_in_rect(rect_num) = bitget(bits(bit), rect_num);
                if any(all(bits_in_rect == 0, 2)) || any(all(bits_in_rect == 1, 2))
                    is_valid = false;
                    break;
                end
            end
        end
    end
    
    if is_valid
        fprintf('找到满足条件的映射！\n');
        disp(bits);
        success = true;
        break;
    end
end

if ~success
    fprintf('未找到满足条件的映射。\n');
end