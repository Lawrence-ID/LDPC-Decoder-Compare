%% 产生LDPC译码模块控制字（使用H矩阵产生），支持NR两种基图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   根据给定的LDPC译码器控制字计算FPGA中ROM的coe文件
%	输入参数：
%       Ix_c ： 校验节点更新控制字，是校验矩阵每个按行块中的每行行重的集合，行向量
%       Iy_c ： 变量节点更新控制字，是校验矩阵每个按列块中每块首列列重的集合，行向量
%       It_a ： 变量节点读取数据的地址交织控制字的集合，行向量
%       It_c ： 变量节点读取数据的数据移位控制字的集合，行向量
%	输出参数：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   clear
   clc
%准备提升因子
%% 0.提升因子Z值表
rom_data = [  2   4   8  16  32  64 128 256 ...
              3   6  12  24  48  96 192 384 ...
              5  10  20  40  80 160 320 ...
              7  14  28  56 112 224 ...
              9  18  36  72 144 288 ...
             11  22  44  88 176 352 ...
             13  26  52 104 208 ...
             15  30  60 120 240];
BG1or2  = 2;                % 基图选择，1->基图1；2->基图2
Ix_c_name_str = strcat('Ix_c_BG',num2str(BG1or2),'_ALL_Zc');
Iyc_Ita_Itc_name_str = strcat('Iyc_Ita_Itc_BG',num2str(BG1or2),'_ALL_Zc');
Ixc_Ixd = [];
buf_Iyc_Ita_Itc = [];
for iii = 1:length(rom_data)
    
    %% 参数给定
    Z       = rom_data(iii);                % 提升因子（并行度）
%     Z       = 6;                % 提升因子（并行度）
    %% 计算LDPC译码器控制字
    % 计算基矩阵（HBG）
    [HBG,rHBG,cHBG] = ldpc_01_HbMatrix_calc(BG1or2,Z);
    % 计算校验矩阵（H）及HBG2特殊位置
    [HMatrix,data_pos,data_spec] = ldpc_02_HMatrix_calc(HBG,Z);
    % 根据码长和码率计算ldpc相关参数
    [LDPC_N, LDPC_K, LDPC_P, LDPC_Q] = ldpc_02_Calc_decparam(BG1or2, Z);
    % 计算LDPC译码器控制字
    [interweave_addr,barrel_shifter,variable_ctrl,check_ctrl] = ldpc_03_Calc_decparam(LDPC_N,LDPC_K,LDPC_P,LDPC_Q,HMatrix, Z);
    % LDPC地址码表列中1的数目的最大值和最小值

    %% 译码器主要参数计算
    % LDPC校验矩阵列中1的数目的最大值和最小值
%     LDPC_COLUMU_MAX_NUM = max(variable_ctrl);%30
%     LDPC_COLUMU_MIN_NUM = min(variable_ctrl);%1
    LDPC_COLUMU_MAX_NUM = 30; %此处应等于BG1的参数值，因为vnu_core和cnu_core是按照最大设计的
    LDPC_COLUMU_MIN_NUM = 1;
    VNU_DLY = LDPC_COLUMU_MAX_NUM + LDPC_COLUMU_MIN_NUM + 4;
    % LDPC校验矩阵行中1的数目的最大值和最小值
%     LDPC_ROW_MAX_NUM = max(check_ctrl);%19
%     LDPC_ROW_MIN_NUM = min(check_ctrl);%3
    LDPC_ROW_MAX_NUM = 19;
    LDPC_ROW_MIN_NUM = 3;
    CNU_DLY = LDPC_ROW_MAX_NUM - LDPC_ROW_MIN_NUM + 10;
    SUM_TIME = length(barrel_shifter) + VNU_DLY + 8;

    %% 处理校验节点更新控制字Ix_c
    %将校验节点更新控制字Ix_c变化为起始数据脉冲和终止数据脉冲，终止脉冲-起始脉冲=校验节点更新控制字
    %起始脉冲Ixc:__|~~|____________|~~|______________|~~|_______________|~~|_____________
    %终止脉冲Ixd:_______________|~~|______________|~~|_______________|~~|________|~~|____
    Ix_c_buf = check_ctrl;
    Ixc = zeros(1,sum(Ix_c_buf));
    Ixd = zeros(1,sum(Ix_c_buf));
    ij = 1;
    for i = 1:length(Ix_c_buf)
        for j = 1:Ix_c_buf(i)
            if j == Ix_c_buf(i)
                Ixc(ij) = 0;
            else
                if j == 1
                    Ixc(ij) = 1;
                else
                    Ixc(ij) = 0;
                end
            end
            ij = ij + 1;
        end
    end
    ij = 1;
    for i = 1:length(Ix_c_buf)
        for j = 1:Ix_c_buf(i)
            if j == Ix_c_buf(i)
                Ixd(ij) = 1;
            else
                Ixd(ij) = 0;
            end
            ij = ij + 1;
        end
    end
    len_Ixc = length(Ixc);
    len_Ixd = length(Ixd);

    Ix_rdvld  = [  ones(1,len_Ixc)    , zeros(1,SUM_TIME - len_Ixc)  ];    %node_ram rden
    Ix_rdaddr = [  0 : 1 : len_Ixc-1  , zeros(1,SUM_TIME - len_Ixc)  ];    %node_ram rdaddr

    Ix_dvld   = [  0,0,ones(1,len_Ixc), zeros(1,SUM_TIME - len_Ixc-2)];    %I_node_dvld,delay from rden
    Ixc       = [  0,0,Ixc zeros(1,SUM_TIME - len_Ixc-2)  ];               %I_node_dpos
    Ixd       = [  0,0,Ixd zeros(1,SUM_TIME - len_Ixd-2)  ];               %I_node_dneg

    Ix_wrvld  = [  zeros(1,CNU_DLY+2) , ones(1,len_Ixc), zeros(1,SUM_TIME - len_Ixc - CNU_DLY - 2) ];    %node_ram wren
    Ix_wraddr = [  zeros(1,CNU_DLY+2) , 0:1:len_Ixc - 1, zeros(1,SUM_TIME - len_Ixc - CNU_DLY - 2) ];    %node_ram wraddr

    %% 将控制字写入到coe文件中
    Ixc_Ixd =...
      Ix_rdvld      * 2^(9 + 1 + 1 + 1 + 1 + 9)...
    + Ix_rdaddr     * 2^(9 + 1 + 1 + 1 + 1)...
    + Ix_dvld       * 2^(9 + 1 + 1 + 1)...
    + Ixc           * 2^(9 + 1 + 1)...
    + Ixd           * 2^(9 + 1)...
    + Ix_wrvld      * 2^(9)...
    + Ix_wraddr;
%     file_Ix_c = strcat(Ix_c_name_str,'.coe');
%     ptr_Ix_c  = fopen(file_Ix_c,'w');
%         fprintf(ptr_Ix_c,'memory_initialization_radix=%d;\n',10);
%         fprintf(ptr_Ix_c,'memory_initialization_vector= \n');
%         fprintf(ptr_Ix_c,'%d,\n',Ixc_Ixd(1:end-1));
%         fprintf(ptr_Ix_c,'%d;\n',Ixc_Ixd(end));
%     fclose(ptr_Ix_c);
    %% 将不同Zc的控制字拼接

    %% 处理变量节点更新控制字Iy_c
    %将变量节点更新控制字Iy_c变化为起始数据脉冲和终止数据脉冲，终止脉冲-起始脉冲=变量节点更新控制字
    %起始脉冲Iyc:__|~~|____________|~~|______________|~~|_______________|~~|_____________
    %终止脉冲Iyd:_______________|~~|______________|~~|_______________|~~|________|~~|____
    Iy_c_buf = variable_ctrl;
    Iyc = zeros(1,sum(Iy_c_buf));
    ij = 1;
    for i = 1:length(Iy_c_buf)
        for j = 1:Iy_c_buf(i)
            if j == Iy_c_buf(i) && Iy_c_buf(i) == 1
                Iyc(ij) = 1;
            else
                if j == Iy_c_buf(i)
                    Iyc(ij) = 0;
                else
                    if j == 1
                        Iyc(ij) = 1;
                    else
                        Iyc(ij) = 0;
                    end
                end
            end
            ij = ij + 1;
        end
    end
    Iyd = circshift(Iyc,-1);

    %% 处理变量节点缓存地址交织控制字
    It_a = interweave_addr;     %地址交织控制字
    %% 处理变量节点缓存数据桶形移位器控制字
    It_c = barrel_shifter;      %正桶型移位控制字
    It_c_inv = It_c;            %反桶型移位控制字
    %% 处理变量节点更新读取LLR控制字
    Iz_c = Iyc;%读取LLR缓存RAM使能
    Iz_a1 = zeros(1,length(Iyc));%读取LLR缓存RAM地址
    ij = 0;
    for i = 1 : length(Iyc)
        if Iyc(i) == 1
            Iz_a1(i) = ij;
            ij = ij + 1;
        else
            Iz_a1(i) = ij - 1;
        end
    end
    len_Iyc = length(Iyc);
    %% 拼接变量节点更新控制字Iy_c、地址正交织控制字It_a、地址反交织控制字It_a_dly、正桶形移位控制字It_c1/It_c3/It_c3、反桶形移位控制字It_c1_dly/It_c3_dly/It_c3_dly
    node_rden   = [  ones(1,len_Iyc),  zeros(1,SUM_TIME - len_Iyc)  ];              %节点缓存读使能
    node_rdaddr = [  It_a           ,  zeros(1,SUM_TIME - len_Iyc)  ];              %节点缓存读地址

    barrel_c1   = [zeros(1,2),It_c,zeros(1,SUM_TIME - len_Iyc - 2)];                %正桶形移位控制字（一级）

    llr_rden    = [zeros(1,4),Iz_c ,zeros(1,SUM_TIME - len_Iyc - 4)];               %LLR缓存读使能
    llr_rdaddr  = [zeros(1,4),Iz_a1,zeros(1,SUM_TIME - len_Iyc - 4)];               %LLR缓存读地址

    llr_dvld    = [zeros(1,5),Iz_c,           zeros(1,SUM_TIME - len_Iyc - 5)];     %I_llr_dvld
    node_dvld   = [zeros(1,5),ones(1,len_Iyc),zeros(1,SUM_TIME - len_Iyc - 5)];     %I_node_dvld
    node_dneg   = [zeros(1,5),Iyd,            zeros(1,SUM_TIME - len_Iyc - 5)];     %I_node_dneg

    barrel_c1_inv = [zeros(1,5+VNU_DLY+0),It_c_inv,       zeros(1,SUM_TIME - len_Iyc - (5+VNU_DLY+0))]; %反桶形移位控制字（一级）
    node_wren     = [zeros(1,5+VNU_DLY+3),ones(1,len_Iyc),zeros(1,SUM_TIME - len_Iyc - (5+VNU_DLY+3))]; %节点缓存写使能
    node_wraddr   = [zeros(1,5+VNU_DLY+3),It_a,           zeros(1,SUM_TIME - len_Iyc - (5+VNU_DLY+3))]; %节点缓存写地址


    %% 拼接变量节点更新控制字Iy_c、地址正交织控制字It_a、地址反交织控制字It_a_dly、正桶形移位控制字It_c1/It_c3/It_c3、反桶形移位控制字It_c1_dly/It_c3_dly/It_c3_dly
    buf_Iyc_Ita_Itc = [buf_Iyc_Ita_Itc;...
      node_rden     * 2^(9 + 1 + 9 + 1 + 1 + 1 + 8 + 1 + 9  + 9)...
    + node_rdaddr   * 2^(9 + 1 + 9 + 1 + 1 + 1 + 8 + 1 + 9)...
    + barrel_c1     * 2^(9 + 1 + 9 + 1 + 1 + 1 + 8 + 1)...
    + llr_rden      * 2^(9 + 1 + 9 + 1 + 1 + 1 + 8)...
    + llr_rdaddr    * 2^(9 + 1 + 9 + 1 + 1 + 1)...
    + llr_dvld      * 2^(9 + 1 + 9 + 1 + 1)...
    + node_dvld     * 2^(9 + 1 + 9 + 1)...
    + node_dneg     * 2^(9 + 1 + 9)...
    + barrel_c1_inv * 2^(9 + 1)...
    + node_wren     * 2^(9)...
    + node_wraddr];
    %% 将不同Zc的控制字拼接

end
    %% 将控制字存储为coe形式
    file_Ix_c = strcat(Ix_c_name_str,'.coe');
    ptr_Ix_c  = fopen(file_Ix_c,'w');
        fprintf(ptr_Ix_c,'memory_initialization_radix=%d;\n',10);
        fprintf(ptr_Ix_c,'memory_initialization_vector= \n');
        fprintf(ptr_Ix_c,'%d,\n',Ixc_Ixd(1:end-1));
        fprintf(ptr_Ix_c,'%d;\n',Ixc_Ixd(end));
    fclose(ptr_Ix_c);
    %% 将控制字存储为coe形式
    file_Iyc_Ita_Itc = strcat(Iyc_Ita_Itc_name_str,'.coe');
    buf_Iyc_Ita_Itc = buf_Iyc_Ita_Itc';
    ptr_Iyc_Ita_Itc  = fopen(file_Iyc_Ita_Itc,'w');
        fprintf(ptr_Iyc_Ita_Itc,'memory_initialization_radix=%d;\n',10);
        fprintf(ptr_Iyc_Ita_Itc,'memory_initialization_vector= \n');
        fprintf(ptr_Iyc_Ita_Itc,'%d,\n',buf_Iyc_Ita_Itc(1:end-1));
        fprintf(ptr_Iyc_Ita_Itc,'%d;\n',buf_Iyc_Ita_Itc(end));
    fclose(ptr_Iyc_Ita_Itc);
    display('COE generate finish!')
