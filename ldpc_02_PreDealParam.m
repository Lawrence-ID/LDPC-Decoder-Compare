%% 编译码参数预处理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入参数
%   srcData  ：信源数据
%   BG1or2   ：选择基图，1：BG1,2：BG2
%   Z        ：提升因子（并行度）
% 输出参数
%   HBG      ：校验基矩阵
%   K        ：信息位长度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [crcLen,K] = ldpc_02_PreDealParam(srcData,BG1or2,Zc)
%% 1.信源输入处理
    if iscolumn(srcData)
        
    else
        srcData = srcData';
    end
    srcDataLen = length(srcData);
    
    
%% 2.计算基图编码块大小（Kcb）
    if BG1or2 == 1
        Kcb = 8448;
    elseif BG1or2 == 2
        Kcb = 3840;
    else
        error('Error!!! \nInput BG1or2(d%) is error!!! \nPlease correct!!!', BG1or2);
    end
    
%% 3.计算编码块数目（cbNum）
    if srcDataLen <= Kcb
        crcLen  = 0;
        cbNum   = 1;
        dataLen = srcDataLen;
    else
        crcLen  = 24;
        cbNum   = ceil( srcDataLen / (Kcb - crcLen) );
        dataLen = srcDataLen + cbNum * crcLen;
    end
    K1 = dataLen / cbNum;
    
%% 2.计算基图编码块大小（Kcb）
    if BG1or2 == 1
        Kb = 22;
    elseif BG1or2 == 2
        if srcDataLen > 640
            Kb = 10;
        elseif srcDataLen > 560
            Kb = 9;
        elseif srcDataLen > 192
            Kb = 8;
        else
            Kb = 6;
        end
    else
        error('Error!!! \nInput BG1or2(d%) is error!!! \nPlease correct!!!', BG1or2);
    end
    
    if BG1or2 == 1
        if Kb * Zc >= K1
            K = 22 * Zc;
        end
    elseif BG1or2 == 2
        if Kb * Zc >= K1
            K = 10 * Zc;
        end
    else
        error('Error!!! \nInput BG1or2(d%) is error!!! \nPlease correct!!!', BG1or2);
    end
    
    
    
    
    
end