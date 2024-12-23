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
Z       = 384;                % 提升因子（并行度）
if BG1or2 == 1
    Rate = 22/68; % 码率
    kb   = 22;
else
    Rate = 10/52;
    kb   = 10;
end
intera  = 5;               % 迭代次数


%% 计算基矩阵（HBG）
[HBG,rHBG,cHBG] = ldpc_01_HbMatrix_calc(BG1or2,Z);

%% 计算校验矩阵（H）及HBG2特殊位置，LDPC译码器控制字
[HMatrix,data_pos,data_spec] = ldpc_02_HMatrix_calc(HBG,Z);
[LDPC_N, LDPC_K, LDPC_P, LDPC_Q] = ldpc_02_Calc_decparam(BG1or2, Z);
[interweave_addr,barrel_shifter,variable_ctrl,check_ctrl, V] = ldpc_03_Calc_decparam(LDPC_N,LDPC_K,LDPC_P,LDPC_Q,HMatrix, Z);

% 打开文件
llrIn_file_lineNum = 384 * 68;
test_input_file_lineNum = 384 * 22;
llrIn_file = fopen("C:\Users\Pngxi\Desktop\llrIn.txt", 'r');
test_input_file = fopen("C:\Users\Pngxi\Desktop\test_input_bit.txt", "r");
    
% 初始化一个空向量存储数据
llr_data = zeros(1, llrIn_file_lineNum);
src_data = zeros(1, test_input_file_lineNum);
    
% 逐行读取前 x 行
for i = 1:llrIn_file_lineNum
    line = fgetl(llrIn_file);  % 读取一行
    if ischar(line)
        llr_data(i) = str2double(line);  % 将字符串转换为数字并存储到数组
    else
        % 如果文件行数少于 x 行，提前结束
        warning('文件行数少于指定的行数 x');
        llr_data = llr_data(1:i-1);  % 截断数据
        break;
    end
end

for i = 1:test_input_file_lineNum
    line = fgetl(test_input_file);  % 读取一行
    if ischar(line)
        src_data(i) = str2double(line);  % 将字符串转换为数字并存储到数组
    else
        % 如果文件行数少于 x 行，提前结束
        warning('文件行数少于指定的行数 x');
        src_data = src_data(1:i-1);  % 截断数据
        break;
    end
end
    
% 关闭文件
fclose(llrIn_file);
fclose(test_input_file);

[dec_data_flooding, llr_flooding] = ldpc_decoder_multrate_serial(llr_data,intera,LDPC_N,LDPC_K, LDPC_P, LDPC_Q,interweave_addr+1,barrel_shifter,variable_ctrl,check_ctrl, Z);
err_bit_flooding = sum(xor(src_data,dec_data_flooding));%一帧错误bit

[dec_data_layered1, llr_layered1]  = ldpc_decoder_layered(llr_data,intera,HBG, HMatrix, check_ctrl, Z);
err_bit_layered1 = sum(xor(src_data,dec_data_layered1));%一帧错误bit

[dec_data_layered2, llr_layered2]  = ldpc_decoder_layered_02(llr_data,intera,HBG, HMatrix, check_ctrl, Z);
err_bit_layered2 = sum(xor(src_data,dec_data_layered2));%一帧错误bit

dec_data_new = nr_ldpc_decoder_new(llr_data,HBG,Z,1,Rate,BG1or2);
err_bit_layered_new = sum(xor(src_data,dec_data_new));%一帧错误bit
