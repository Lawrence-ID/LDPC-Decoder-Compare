%% 5G NR LDPC������������ ������ROM_A���ɴ���
% BG1:ROM_A:4*22
% BG2:ROM_A:4*10
% ldpc_01_gen_shift_romA_coe
% HBGԪ����Ϊ-1��ȫ�㷽����Ϊ0������λ����Ϊ�����������λ��ѭ�����ƴ�ֵ��
% ���д洢��4·����
clear
clc
close all


%% 0.��������Zֵ��
rom_data = [  2   4   8  16  32  64 128 256 ...
              3   6  12  24  48  96 192 384 ...
              5  10  20  40  80 160 320 ...
              7  14  28  56 112 224 ...
              9  18  36  72 144 288 ...
             11  22  44  88 176 352 ...
             13  26  52 104 208 ...
             15  30  60 120 240];
         

%% 1.����BG1����Zֵ�ľ���Aѭ����λϵ��
temp_bg1 = [];
BG1or2 = 1; % ��ͼѡ��1->��ͼ1��2->��ͼ2
for ii = 1:51
    Z = rom_data(ii);
    %���������HBG��
    [HBG,rHBG,cHBG] = ldpc_01_HbMatrix_calc(BG1or2,Z);
    HBG = HBG + 1;
    temp_bg1 = [temp_bg1 HBG(1:4,1:22)];
end
temp_bg1 = temp_bg1';
romA_shiftsize_BG1 =...
      temp_bg1(:,1)	.* 2^(9 + 9 + 9)...
    + temp_bg1(:,2)	.* 2^(9 + 9)...
    + temp_bg1(:,3)	.* 2^(9)...
    + temp_bg1(:,4);


%% 2.����BG2����Zֵ�ľ���Aѭ����λϵ��
temp_bg2 = [];
BG1or2 = 2; % ��ͼѡ��1->��ͼ1��2->��ͼ2
for ii = 1:51
    Z = rom_data(ii);
    %���������HBG��
    [HBG,rHBG,cHBG] = ldpc_01_HbMatrix_calc(BG1or2,Z);
    HBG = HBG + 1;
    temp_bg2 = [temp_bg2 HBG(1:4,1:10) zeros(4,12)];
end
temp_bg2 = temp_bg2';
romA_shiftsize_BG2 =...
      temp_bg2(:,1)	.* 2^(9 + 9 + 9)...
    + temp_bg2(:,2)	.* 2^(9 + 9)...
    + temp_bg2(:,3)	.* 2^(9)...
    + temp_bg2(:,4);


%% ��ѭ����λϵ��д�뵽COE��
write_data = [romA_shiftsize_BG1;romA_shiftsize_BG2];
%׼����������
coe_name_str = 'ldpc_enc_BG_romA';
%����������д�뵽COE��
fpga_generate_coe(write_data,coe_name_str)


