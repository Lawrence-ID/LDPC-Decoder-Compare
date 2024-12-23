
%% ����У�������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������
%   HBG  ��������
%   Z    ���������ӣ����жȣ�
%   BG1or2 ����ͼѡ��1->��ͼ1��2->��ͼ2
% �������
%   HMatrix ��У�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HMatrix,data_pos,data_spec] = ldpc_02_HMatrix_calc(HBG,Z)
    
    %% 1.��չ������HBG��ΪУ�����H��
    [r_Hb,c_Hb] = size(HBG);
    data_pos = find(HBG(:,c_Hb - r_Hb+1) ~= -1);
    data_spec = HBG(data_pos,c_Hb - r_Hb+1);
    
    HMatrix = zeros(r_Hb*Z,c_Hb*Z);
    for ii = 1 : r_Hb
        for jj = 1 : c_Hb
            if HBG(ii,jj) == -1
                HMatrix((ii-1)*Z+1:(ii-1)*Z+Z,(jj-1)*Z+1:(jj-1)*Z+Z) = zeros(Z,Z);
            elseif HBG(ii,jj) == 0
                HMatrix((ii-1)*Z+1:(ii-1)*Z+Z,(jj-1)*Z+1:(jj-1)*Z+Z) = eye(Z);
            else
                HMatrix((ii-1)*Z+1:(ii-1)*Z+Z,(jj-1)*Z+1:(jj-1)*Z+Z) = circshift(eye(Z),[0,HBG(ii,jj)]);
            end
        end
    end
end