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
intera  = 8;               % 迭代次数
llr_width = 6;             % llr位宽
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

bar = waitbar(0, "Processing...");

for idx = start_idx : end_idx
    for frame = 1 : frame_num(idx)
        progress = ((idx - 1) * frame_num(idx) + frame) / (length(EbN0) * frame_num(idx));
        waitbar(progress, bar, sprintf('Progress: idx: %d/%d, frame: %d/%d', idx, end_idx - start_idx + 1, frame, frame_num(idx)));
        %% 产生1010信源
            Len_1010 = (cHBG - rHBG) * Z; % 信源比特长度 = (#cols(HBG)-#rows(HBG))*Z
            %srcData = reshape([ones(1,Len_1010/2);zeros(1,Len_1010/2)],1,Len_1010);% 1010序列
            srcData = randi([0 1 ],1,Len_1010); % 随机信源

        %% 进行LDPC编码
            % LDPC编码方法一
            % encData = ldpc_04_Encoder(HBG,Z,srcData); % 码字比特长度 = #cols(HBG)*Z

            % LDPC编码方法二：%近似下三角矩阵编码的方法，求B_inv
            % load('srcData_crc24.mat'); srcData = srcData(1:kb*Z);
            [encoded_bits, H, Z_c, encoded_bits_original] = ldpc_encode(srcData, BG1or2,kb);

        %% 将编码数据写到TXT，与FPGA对比验证
    %             aa = reshape(encData',Z,[]);
    %             aa(4,1:end) = 0;
    %             bb = reshape(aa,1,[]);
    %             cc = reshape(bb,4,[]);
    %             dd = bi2de(cc','left-msb');
    %             ee = dec2hex(dd);
    %             ff = reshape(ee',4,[]);
    %             gg = ff';
    %             %写编码输出数据到txt
    %             coe_name_str = strcat('bg',num2str(BG1or2-1),'_','zc',num2str(Z),'_in')
    %             fpga_generate_txt(encData',coe_name_str);
    %             %写编码输出数据到txt
    %             coe_name_str = strcat('bg',num2str(BG1or2-1),'_','zc',num2str(Z),'_in_kb',num2str(kb))
    %             fpga_generate_txt(encoded_bits_original,coe_name_str);
    %         end
    % 
    %         plot(bi2de(reshape(encData(length(srcData)+1:end),4,[])','left-msb'))
    %         a = (bi2de(reshape(encData(length(srcData)+1:end),4,[])','left-msb'))';
    %         a(1:Z/4:end)


%     encoded_bits_original(1:Z*2,1) = 0;
    
    
       %% 信道映射及加噪
            %txData = 1 - 2 * encData; %BPSK bit to symbol conversion
            txData = 1 - 2 * encoded_bits_original';
            en = 10^(EbN0(idx)/10);
            sigma = sqrt(1/(2*Rate*en));  
            rxData = txData + sigma * randn(1,length(txData)); %AWGN channel I
            SNR_dB(idx) = snr(txData, rxData - txData);
            %rxData = txData;
       %% 接收及定点
            rxData(rxData > 4) = 4;
            rxData(rxData < -4) = -4;
            llr_data = floor(rxData*2^4);

       %% LDPC全码率并行译码W
            [dec_data_flooding, llr_flooding] = ldpc_decoder_multrate_serial(llr_data,llr_width,intera,LDPC_N,LDPC_K, LDPC_P, LDPC_Q,interweave_addr+1,barrel_shifter,variable_ctrl,check_ctrl, Z);
            err_bit_flooding = sum(xor(srcData,dec_data_flooding));%一帧错误bit
       %% LDCP layered decoding
            [dec_data_layered, llr_layered, it]  = ldpc_decoder_layered_02(llr_data, llr_width, intera, HBG, HMatrix, check_ctrl, Z);
            err_bit_layered = sum(xor(srcData,dec_data_layered));%一帧错误bit
       %% 错误统计
            errs_flooding(idx) = errs_flooding(idx) + err_bit_flooding;%总错误bit
            [frame errs_flooding(idx)/(LDPC_K*frame)]; %误码率

            errs_layered(idx) = errs_layered(idx) + err_bit_layered;%总错误bit
            [frame errs_layered(idx)/(LDPC_K*frame)]; %误码率
    end
    %% Statistic
    BER_flooding(idx) = errs_flooding(idx)/(LDPC_K*frame_num(idx)); % Bit Error Rate
    BER_layered(idx) = errs_layered(idx)/(LDPC_K*frame_num(idx));
end

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
filename = sprintf('result_BG%d_Z%d_%dframes_%diteras_fix241223.txt', BG1or2, Z, max(frame_num(start_idx:end_idx)), intera);
fileID = fopen(filename, 'w');

% 写入数据
fprintf(fileID, 'Eb/N0\tBER_flooding\tBER_layered\n'); % 写入标题
for i = 1:length(EbN0)
    fprintf(fileID, '%f\t%e\t%e\n', EbN0(i), BER_flooding(i), BER_layered(i)); % 写入数据
end

% 关闭文件
fclose(fileID);


