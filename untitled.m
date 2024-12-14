% 参数设置
M = 64; % 16QAM
num_trials = 1e6; % 随机尝试次数
bit_num = log2(M);
% 星座点排列（以坐标编号为准）
rectangles = [1, 2, 5, 6; 2, 3, 6, 7; 3, 4, 7, 8; 
              5, 6, 9, 10; 6, 7, 10, 11; 7, 8, 11, 12;
              9, 10, 13, 14; 10, 11, 14, 15; 11, 12, 15, 16];

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
        fprintf('找到满足条件的映射！\n');
        disp(bits);
        success = true;
        break;
    end
end

if ~success
    fprintf('未找到满足条件的映射。\n');
end