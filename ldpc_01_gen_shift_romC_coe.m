%% 5G NR LDPC编码提升因子 控制字ROM_A生成代码
% BG1:ROM_C:42*22
% BG2:ROM_C:38*10
% ldpc_01_gen_shift_romC_coe
% HBG元素若为-1：全零方阵，若为0，代表单位阵，若为正数，则代表单位阵循环右移此值次
% 按列存储，42路并行
clear
clc
close all


%% 0.提升因子Z值表
rom_data = [  2   4   8  16  32  64 128 256 ...
              3   6  12  24  48  96 192 384 ...
              5  10  20  40  80 160 320 ...
              7  14  28  56 112 224 ...
              9  18  36  72 144 288 ...
             11  22  44  88 176 352 ...
             13  26  52 104 208 ...
             15  30  60 120 240];
         
         
%% 1.计算BG1所有Z值的矩阵A循环移位系数
temp_bg1 = [];
BG1or2 = 1; % 基图选择，1->基图1；2->基图2
for ii = 1:51
    Z = rom_data(ii);
    %计算基矩阵（HBG）
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

%% 2.计算BG2所有Z值的矩阵A循环移位系数
temp_bg2 = [];
BG1or2 = 2; % 基图选择，1->基图1；2->基图2
for ii = 1:51
    Z = rom_data(ii);
    %计算基矩阵（HBG）
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

%% 将循环移位系数写入到COE中
write_data = [romC_shiftsize_BG1;romC_shiftsize_BG2];
%准备提升因子
coe_name_str = 'ldpc_enc_BG_romC';
%将提升因子写入到COE中
fpga_generate_coe(write_data,coe_name_str)
% fpga_generate_coe_binary(write_data,coe_name_str)


