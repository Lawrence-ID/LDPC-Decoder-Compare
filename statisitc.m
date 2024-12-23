BG1or2  = 1;                % 基图选择，1->基图1；2->基图2
Z       = 384;               % 提升因子（并行度）
if BG1or2 == 1
    Rate = 22/68; % 码率
    kb   = 22;
else
    Rate = 10/52;
    kb   = 10;
end
Zc = 384;
[LDPC_N, LDPC_K, LDPC_P, LDPC_Q] = ldpc_02_Calc_decparam(BG1or2, Z);

[intera10, EbN0, BER_flooding_10iteras, BER_layered_10iteras, FER_flooding_10iteras, FER_layered_10iteras] = read_statistic('result_BG1_Z384_10000frames_10iteras.txt', LDPC_K);
[intera8, EbN0, BER_flooding_8iteras, BER_layered_8iteras, FER_flooding_8iteras, FER_layered_8iteras]      = read_statistic('result_BG1_Z384_10000frames_8iteras.txt', LDPC_K);
[intera5, EbN0, BER_flooding_5iteras, BER_layered_5iteras, FER_flooding_5iteras, FER_layered_5iteras]      = read_statistic('result_BG1_Z384_10000frames_5iteras.txt', LDPC_K);
[intera6, ~, BER_flooding_6iteras, BER_layered_6iteras, FER_flooding_6iteras, FER_layered_6iteras]         = read_statistic('result_BG1_Z384_200frames_6iteras.txt', LDPC_K);

% 创建新的图形窗口
figure;

% 第一个子图：绘制 BER
subplot(1, 2, 1); % 1 行 2 列，第一个子图
semilogy(EbN0, BER_flooding_8iteras, 'b-*', EbN0, BER_layered_8iteras, 'b-o', EbN0, BER_flooding_5iteras, 'r-*', EbN0, BER_layered_5iteras, 'r-o');
xlabel('Eb/N0 (dB)');
ylabel('BER');
title(['Bit Error Rate']);
legend('Flooding(8iteras)', 'Layered(8iteras)', 'Flooding(5iteras)', 'Layered(5iteras)');  

% 第二个子图：绘制 FER
subplot(1, 2, 2); % 1 行 2 列，第二个子图
semilogy(EbN0, FER_flooding_8iteras, 'b-*', EbN0, FER_layered_8iteras, 'b-o', EbN0, FER_flooding_5iteras, 'r-*', EbN0, FER_layered_5iteras, 'r-o');
xlabel('Eb/N0 (dB)');
ylabel('FER');
title(['Frame Error Rate']);
legend('Flooding(8iteras)', 'Layered(8iteras)', 'Flooding(5iteras)', 'Layered(5iteras)');
