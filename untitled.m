figure;
semilogy(SNR_dB,error_rate_MMSE,'b-s',SNR_dB,error_rate_box,'g-x',SNR_dB,error_rate_KB,'r-o');
legend('MMSE','box,k=4','KB,k=4');
title('BER of  MMSE, box and KB');
xlabel('SNR in dB');
ylabel('BER');

figure; % 创建新图形窗口
plot(SNR_dB, PED_count_box, '-o'); % 绘制第一组数据
hold on; % 保持当前图形
plot(SNR_dB, PED_count_KB, '-x'); % 绘制第二组数据

% 添加图例
legend('box,k=4','KB,k=4');

% 添加标题和轴标签
title('PED times on different SNR');
xlabel('SNR (dB)');
ylabel('PED calculation tiumes');

% 显示网格
grid on;