%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入参数
%   BG1or2 ：基图选择，1->基图1；2->基图2
%   Zc  ：提升因子
% 输出参数
% LDPC_N    : 码长
% LDPC_K    : 信息位长度
% LDPC_P    : 校验位长度
% LDPC_Q    : q值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 根据码长和码率计算ldpc相关参数
function [LDPC_N, LDPC_K, LDPC_P, LDPC_Q] = ldpc_02_Calc_decparam(BG1or2, Z_c)
    if(BG1or2 == 1)             %如果基图选择1
        LDPC_K = 22 * Z_c;
        LDPC_N = 68 * Z_c;
    elseif (BG1or2 == 2)        %如果基图选择2
        LDPC_K = 10 * Z_c;
        LDPC_N = 52 * Z_c;
    else
      error('Wrong base graph index in ldpc encoding.');
    end
    LDPC_P = LDPC_N - LDPC_K;
    LDPC_Q = (LDPC_N - LDPC_K) / Z_c;
%     LDPC_Q = Z_c;
end