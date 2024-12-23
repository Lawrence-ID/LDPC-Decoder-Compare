%% 5G NR LDPC Simulation
clear
clc
close all
%׼����������
%% 0.��������Zֵ��
% rom_data = [  2   4   8  16  32  64 128 256 ...
%               3   6  12  24  48  96 192 384 ...
%               5  10  20  40  80 160 320 ...
%               7  14  28  56 112 224 ...
%               9  18  36  72 144 288 ...
%              11  22  44  88 176 352 ...
%              13  26  52 104 208 ...
%              15  30  60 120 240];
         
%% ��������
% for i = 1:length(rom_data)
BG1or2  = 1;                % ��ͼѡ��1->��ͼ1��2->��ͼ2
Z       = 3;                % �������ӣ����жȣ�
if BG1or2 == 1
    Rate = 22/68; % ����
    kb   = 22;
else
    Rate = 10/52;
    kb   = 10;
end
intera  = 8;               % ��������
% frame_num   = [1e3    1e3       1e3     1e3     1e5];
% idx           1    2   3   4   5   6  7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
EbN0        = [0.5 0.6 0.8 0.9 1.0 1.2 1.4 1.5 1.6 1.7 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8 4.0 4.2 4.4 4.6 4.8 5.0];
frame_num   = [1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2 2e2 2e2 2e2 2e2 5e2 5e2 1e3 1e3 2e3 3e3 3e3 5e3 5e3 1e4];
% frame_num   = ones(size(EbN0)) * 1000;
SNR_dB      = zeros(size(EbN0));
start_idx   = 1;
end_idx     = 17;

%% ���������HBG��
[HBG,rHBG,cHBG] = ldpc_01_HbMatrix_calc(BG1or2,Z);

%% ����У�����H����HBG2����λ�ã�LDPC������������
[HMatrix,data_pos,data_spec] = ldpc_02_HMatrix_calc(HBG,Z);
[LDPC_N, LDPC_K, LDPC_P, LDPC_Q] = ldpc_02_Calc_decparam(BG1or2, Z);
[interweave_addr,barrel_shifter,variable_ctrl,check_ctrl, V] = ldpc_03_Calc_decparam(LDPC_N,LDPC_K,LDPC_P,LDPC_Q,HMatrix, Z);

start_time = clock;

errs_flooding = zeros(1,length(EbN0));
BER_flooding = zeros(1,length(EbN0));

errs_layered = zeros(1,length(EbN0));
BER_layered = zeros(1,length(EbN0));


llr_data = mod(0:cHBG * Z - 1, 64);
% ��llr_dataת��Ϊ6λ�з�����
signed_llr_data = llr_data;
signed_llr_data(llr_data > 31) = llr_data(llr_data > 31) - 64

%% LDCP layered decoding
[dec_data_layered, llr_layered, it]  = ldpc_decoder_layered_02(signed_llr_data, intera, HBG, HMatrix, check_ctrl, Z);
err_bit_layered = sum(xor(srcData,dec_data_layered));%һ֡����bit
%% ����ͳ��
errs_flooding(idx) = errs_flooding(idx) + err_bit_flooding;%�ܴ���bit
[frame errs_flooding(idx)/(LDPC_K*frame)]; %������
errs_layered(idx) = errs_layered(idx) + err_bit_layered;%�ܴ���bit
[frame errs_layered(idx)/(LDPC_K*frame)]; %������
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

% ���ļ�����д��
filename = sprintf('result_BG%d_Z%d_%dframes_%diteras.txt', BG1or2, Z, max(frame_num(start_idx:end_idx)), intera);
fileID = fopen(filename, 'w');

% д������
fprintf(fileID, 'Eb/N0\tBER_flooding\tBER_layered\n'); % д�����
for i = 1:length(EbN0)
    fprintf(fileID, '%f\t%e\t%e\n', EbN0(i), BER_flooding(i), BER_layered(i)); % д������
end

% �ر��ļ�
fclose(fileID);


