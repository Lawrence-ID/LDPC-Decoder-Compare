%% 5G NR LDPC������������ ������ROM_A���ɴ���
% BG1:ROM_C:42*22
% BG2:ROM_C:38*10
% ldpc_01_gen_shift_romC_coe
% HBGԪ����Ϊ-1��ȫ�㷽����Ϊ0������λ����Ϊ�����������λ��ѭ�����ƴ�ֵ��
% ���д洢��42·����
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
    temp_bg1 = [temp_bg1 HBG(5:46,1:22)];
end
temp_bg1 = temp_bg1';

[r,c] = size(temp_bg1);
data0 = zeros(r*c,9);
for ii = 1:r
    for jj = 1:c
        data0( ( ii-1)*c + jj,:) = de2bi(temp_bg1(ii,jj),9,'left-msb');
    end
end
data1 = data0';
data2 = reshape(data1,42*9,1122);
romC_shiftsize_BG1 = data2';

%% 2.����BG2����Zֵ�ľ���Aѭ����λϵ��
temp_bg2 = [];
BG1or2 = 2; % ��ͼѡ��1->��ͼ1��2->��ͼ2
for ii = 1:51
    Z = rom_data(ii);
    %���������HBG��
    [HBG,rHBG,cHBG] = ldpc_01_HbMatrix_calc(BG1or2,Z);
    HBG = HBG + 1;
    temp_bg2 = [temp_bg2 [ [HBG(5:42,1:10) zeros(38,12)]; zeros(4,22)] ];
end
temp_bg2 = temp_bg2';

[r,c] = size(temp_bg2);
data0 = zeros(r*c,9);
for ii = 1:r
    for jj = 1:c
        data0( ( ii-1)*c + jj,:) = de2bi(temp_bg2(ii,jj),9,'left-msb');
    end
end
data1 = data0';
data2 = reshape(data1,42*9,1122);
romC_shiftsize_BG2 = data2';

%% ��ѭ����λϵ��д�뵽COE��
write_data = [romC_shiftsize_BG1;romC_shiftsize_BG2];
%׼����������
coe_name_str = 'ldpc_enc_BG_romC';
%����������д�뵽COE��
fpga_generate_coe(write_data,coe_name_str)
% fpga_generate_coe_binary(write_data,coe_name_str)


