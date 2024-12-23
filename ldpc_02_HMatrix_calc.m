
%% 计算校验基矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入参数
%   HBG  ：基矩阵
%   Z    ：提升因子（并行度）
%   BG1or2 ：基图选择，1->基图1；2->基图2
% 输出参数
%   HMatrix ：校验矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HMatrix,data_pos,data_spec] = ldpc_02_HMatrix_calc(HBG,Z)
    
    %% 1.扩展基矩阵（HBG）为校验矩阵（H）
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