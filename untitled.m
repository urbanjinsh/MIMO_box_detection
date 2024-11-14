figure; % 创建新图形窗口
plot(SNR_dB, PED_count_SD, '-o'); % 绘制第一组数据
hold on; % 保持当前图形
plot(SNR_dB, PED_count_ML, '-x'); % 绘制第二组数据

% 添加图例
legend('SD','ML');

% 添加标题和轴标签
title('PED times on different SNR');
xlabel('SNR (dB)');
ylabel('PED calculation tiumes');

% 显示网格
grid on;

figure;
semilogy(SNR_dB,error_rate_MMSE,'b-s',SNR_dB,error_rate_SD,'g-x',SNR_dB,error_rate_ML,'r-o');
legend('MMSE','SD','ML');
title('BER of  MMSE, SD and ML');
xlabel('SNR in dB');
ylabel('BER')