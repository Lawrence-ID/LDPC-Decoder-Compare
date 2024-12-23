%% 5G NR LDPC Simulation
clear
clc
close all
%准备提升因子
%% 0.提升因子Z值表
% rom_data = [  2   4   8  16  32  64 128 256 ...
%               3   6  12  24  48  96 192 384 ...
%               5  10  20  40  80 160 320 ...
%               7  14  28  56 112 224 ...
%               9  18  36  72 144 288 ...
%              11  22  44  88 176 352 ...
%              13  26  52 104 208 ...
%              15  30  60 120 240];
         
%% 给定参数
% for i = 1:length(rom_data)
BG1or2  = 1;                % 基图选择，1->基图1；2->基图2
Z       = 3;                % 提升因子（并行度）
if BG1or2 == 1
    Rate = 22/68; % 码率
    kb   = 22;
else
    Rate = 10/52;
    kb   = 10;
end
intera  = 8;               % 迭代次数
% frame_num   = [1e3    1e3       1e3     1e3     1e5];
% idx           1    2   3   4   5   6  7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
EbN0        = [0.5 0.6 0.8 0.9 1.0 1.2 1.4 1.5 1.6 1.7 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8 4.0 4.2 4.4 4.6 4.8 5.0];
frame_num   = [1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 2e2 2e2 2e2 2e2 5e2 5e2 1e3 1e3 2e3 3e3 3e3 5e3 5e3 1e4];
% frame_num   = ones(size(EbN0)) * 1000;
SNR_dB      = zeros(size(EbN0));
start_idx   = 1;
end_idx     = 17;

%% 计算基矩阵（HBG）
[HBG,rHBG,cHBG] = ldpc_01_HbMatrix_calc(BG1or2,Z);

%% 计算校验矩阵（H）及HBG2特殊位置，LDPC译码器控制字
[HMatrix,data_pos,data_spec] = ldpc_02_HMatrix_calc(HBG,Z);
[LDPC_N, LDPC_K, LDPC_P, LDPC_Q] = ldpc_02_Calc_decparam(BG1or2, Z);
[interweave_addr,barrel_shifter,variable_ctrl,check_ctrl, V] = ldpc_03_Calc_decparam(LDPC_N,LDPC_K,LDPC_P,LDPC_Q,HMatrix, Z);

start_time = clock;

errs_flooding = zeros(1,length(EbN0));
BER_flooding = zeros(1,length(EbN0));

errs_layered = zeros(1,length(EbN0));
BER_layered = zeros(1,length(EbN0));


llr_data = mod(0:cHBG * Z - 1, 64);
% 将llr_data转换为6位有符号数
signed_llr_data = llr_data;
signed_llr_data(llr_data > 31) = llr_data(llr_data > 31) - 64

%% LDCP layered decoding
[dec_data_layered, llr_layered, it]  = ldpc_decoder_layered_02(signed_llr_data, intera, HBG, HMatrix, check_ctrl, Z);
err_bit_layered = sum(xor(srcData,dec_data_layered));%一帧错误bit
%% 错误统计
errs_flooding(idx) = errs_flooding(idx) + err_bit_flooding;%总错误bit
[frame errs_flooding(idx)/(LDPC_K*frame)]; %误码率
errs_layered(idx) = errs_layered(idx) + err_bit_layered;%总错误bit
[frame errs_layered(idx)/(LDPC_K*frame)]; %误码率
%% Statistic
BER_flooding(idx) = errs_flooding(idx)/(LDPC_K*frame_num(idx)); % Bit Error Rate

end_time = clock;
etime(end_time,start_time)

EbN0
SNR_dB
BER_flooding
BER_layered

BER_flooding(idx) = errs_flooding(idx)/(LDPC_K*frame_num(idx));
FER_flooding = 1 - (1 - BER_flooding) .^ LDPC_K;

BER_layered(idx) = errs_layered(idx)/(LDPC_K*frame_num(idx));
FER_layered = 1 - (1 - BER_layered) .^ LDPC_K;

semilogy(EbN0,BER_flooding,'b-*', EbN0, BER_layered,'r-*');
% semilogy(EbN0,FER_flooding,'b-*', EbN0, FER_layered,'r-*');
close(bar);

% 打开文件进行写入
filename = sprintf('result_BG%d_Z%d_%dframes_%diteras.txt', BG1or2, Z, max(frame_num(start_idx:end_idx)), intera);
fileID = fopen(filename, 'w');

% 写入数据
fprintf(fileID, 'Eb/N0\tBER_flooding\tBER_layered\n'); % 写入标题
for i = 1:length(EbN0)
    fprintf(fileID, '%f\t%e\t%e\n', EbN0(i), BER_flooding(i), BER_layered(i)); % 写入数据
end

% 关闭文件
fclose(fileID);


