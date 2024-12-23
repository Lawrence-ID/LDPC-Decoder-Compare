%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������
%   BG1or2 ����ͼѡ��1->��ͼ1��2->��ͼ2
%   Zc  ����������
% �������
% LDPC_N    : �볤
% LDPC_K    : ��Ϣλ����
% LDPC_P    : У��λ����
% LDPC_Q    : qֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �����볤�����ʼ���ldpc��ز���
function [LDPC_N, LDPC_K, LDPC_P, LDPC_Q] = ldpc_02_Calc_decparam(BG1or2, Z_c)
    if(BG1or2 == 1)             %�����ͼѡ��1
        LDPC_K = 22 * Z_c;
        LDPC_N = 68 * Z_c;
    elseif (BG1or2 == 2)        %�����ͼѡ��2
        LDPC_K = 10 * Z_c;
        LDPC_N = 52 * Z_c;
    else
      error('Wrong base graph index in ldpc encoding.');
    end
    LDPC_P = LDPC_N - LDPC_K;
    LDPC_Q = (LDPC_N - LDPC_K) / Z_c;
%     LDPC_Q = Z_c;
end