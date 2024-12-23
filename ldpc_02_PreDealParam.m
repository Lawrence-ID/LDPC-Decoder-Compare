%% ���������Ԥ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������
%   srcData  ����Դ����
%   BG1or2   ��ѡ���ͼ��1��BG1,2��BG2
%   Z        ���������ӣ����жȣ�
% �������
%   HBG      ��У�������
%   K        ����Ϣλ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [crcLen,K] = ldpc_02_PreDealParam(srcData,BG1or2,Zc)
%% 1.��Դ���봦��
    if iscolumn(srcData)
        
    else
        srcData = srcData';
    end
    srcDataLen = length(srcData);
    
    
%% 2.�����ͼ������С��Kcb��
    if BG1or2 == 1
        Kcb = 8448;
    elseif BG1or2 == 2
        Kcb = 3840;
    else
        error('Error!!! \nInput BG1or2(d%) is error!!! \nPlease correct!!!', BG1or2);
    end
    
%% 3.����������Ŀ��cbNum��
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
    
%% 2.�����ͼ������С��Kcb��
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