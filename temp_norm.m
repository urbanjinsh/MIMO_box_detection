clear all;
close all;

% 定义接收信号向量 y
y = [1 + 1i; 2 + 2i]; % 示例接收信号向量

% 定义信道矩阵 H
H = [1, 0; 0, 1]; % 示例信道矩阵

% 定义发送符号向量 temp_x
temp_x = [1 + 1i; -1 - 1i]; % 示例发送符号向量

% 计算接收信号与通过信道矩阵传输的符号之间的距离
temp = norm(y - H * temp_x);

% 显示计算结果
disp(['距离: ', num2str(temp)]);

disp(['距离: ', H(0)]);