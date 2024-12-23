clear
clc
close all

HBG=[1 -1 2 1; 0 2 -1 0;2 -1 1 1]; Zc = 3; kb=1; LDPC_K = kb * Zc;LDPC_Q = 3;LDPC_N = Zc * (kb + LDPC_Q);
[HMatrix, ~, ~] = ldpc_02_HMatrix_calc(HBG, Zc);LDPC_P = Zc * LDPC_Q;
llr_data = [-10,-34,45,16,-34,30,-60,2,-10,20,-17,-39];
intera = 2;

% frame_num   = [1e3    1e3       1e3     1e3     1e5     1e6];
% EbN0        = [1.2    1.4       1.6     1.8     2.2     2.4];
% SNR_dB      = [0 0 0 0 0];

% [rHBG, cHBG] = size(HBG);
% [interweave_addr,barrel_shifter,variable_ctrl,check_ctrl, V] = ldpc_03_Calc_decparam(LDPC_N,LDPC_K,LDPC_P,LDPC_Q,HMatrix, Z);
% start_time = clock;
% errs = zeros(1,length(EbN0));
% ferrs = zeros(1,length(EbN0));
% ber = zeros(1,length(EbN0));
% fer = zeros(1, length(EbN0));

errs_layered = zeros(1,length(EbN0));
ferrs_layered = zeros(1,length(EbN0));
ber_layered = zeros(1,length(EbN0));
fer_layered = zeros(1, length(EbN0));

% Len_1010 = (cHBG - rHBG) * Z; % 信源比特长度 = (#cols(HBG)-#rows(HBG))*Z
% %srcData = reshape([ones(1,Len_1010/2);zeros(1,Len_1010/2)],1,Len_1010);% 1010序列
% srcData = randi([0 1 ],1,Len_1010); % 随机信源
% 
% encData = ldpc_04_Encoder(HBG,Z,srcData); % 码字比特长度 = #cols(HBG)*Z

[dec_data_flooding1, llr_flooding1] = ldpc_decoder_multrate_serial(llr_data,1,LDPC_N,LDPC_K, LDPC_P, LDPC_Q,interweave_addr+1,barrel_shifter,variable_ctrl,check_ctrl, Z);
[dec_data_flooding2, llr_flooding2] = ldpc_decoder_multrate_serial(llr_data,2,LDPC_N,LDPC_K, LDPC_P, LDPC_Q,interweave_addr+1,barrel_shifter,variable_ctrl,check_ctrl, Z);
[dec_data_flooding3, llr_flooding3] = ldpc_decoder_multrate_serial(llr_data,3,LDPC_N,LDPC_K, LDPC_P, LDPC_Q,interweave_addr+1,barrel_shifter,variable_ctrl,check_ctrl, Z);

[dec_data_layered1, llr_layered1]  = ldpc_decoder_layered(llr_data,1,HBG, HMatrix, check_ctrl, Z);
[dec_data_layered2, llr_layered2]  = ldpc_decoder_layered(llr_data,2,HBG, HMatrix, check_ctrl, Z);