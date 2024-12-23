%% ����LDPC����ģ������֣�ʹ��H�����������֧��NR���ֻ�ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ���ݸ�����LDPC�����������ּ���FPGA��ROM��coe�ļ�
%	���������
%       Ix_c �� У��ڵ���¿����֣���У�����ÿ�����п��е�ÿ�����صļ��ϣ�������
%       Iy_c �� �����ڵ���¿����֣���У�����ÿ�����п���ÿ���������صļ��ϣ�������
%       It_a �� �����ڵ��ȡ���ݵĵ�ַ��֯�����ֵļ��ϣ�������
%       It_c �� �����ڵ��ȡ���ݵ�������λ�����ֵļ��ϣ�������
%	���������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   clear
   clc
%׼����������
%% 0.��������Zֵ��
rom_data = [  2   4   8  16  32  64 128 256 ...
              3   6  12  24  48  96 192 384 ...
              5  10  20  40  80 160 320 ...
              7  14  28  56 112 224 ...
              9  18  36  72 144 288 ...
             11  22  44  88 176 352 ...
             13  26  52 104 208 ...
             15  30  60 120 240];
BG1or2  = 2;                % ��ͼѡ��1->��ͼ1��2->��ͼ2
Ix_c_name_str = strcat('Ix_c_BG',num2str(BG1or2),'_ALL_Zc');
Iyc_Ita_Itc_name_str = strcat('Iyc_Ita_Itc_BG',num2str(BG1or2),'_ALL_Zc');
Ixc_Ixd = [];
buf_Iyc_Ita_Itc = [];
for iii = 1:length(rom_data)
    
    %% ��������
    Z       = rom_data(iii);                % �������ӣ����жȣ�
%     Z       = 6;                % �������ӣ����жȣ�
    %% ����LDPC������������
    % ���������HBG��
    [HBG,rHBG,cHBG] = ldpc_01_HbMatrix_calc(BG1or2,Z);
    % ����У�����H����HBG2����λ��
    [HMatrix,data_pos,data_spec] = ldpc_02_HMatrix_calc(HBG,Z);
    % �����볤�����ʼ���ldpc��ز���
    [LDPC_N, LDPC_K, LDPC_P, LDPC_Q] = ldpc_02_Calc_decparam(BG1or2, Z);
    % ����LDPC������������
    [interweave_addr,barrel_shifter,variable_ctrl,check_ctrl] = ldpc_03_Calc_decparam(LDPC_N,LDPC_K,LDPC_P,LDPC_Q,HMatrix, Z);
    % LDPC��ַ�������1����Ŀ�����ֵ����Сֵ

    %% ��������Ҫ��������
    % LDPCУ���������1����Ŀ�����ֵ����Сֵ
%     LDPC_COLUMU_MAX_NUM = max(variable_ctrl);%30
%     LDPC_COLUMU_MIN_NUM = min(variable_ctrl);%1
    LDPC_COLUMU_MAX_NUM = 30; %�˴�Ӧ����BG1�Ĳ���ֵ����Ϊvnu_core��cnu_core�ǰ��������Ƶ�
    LDPC_COLUMU_MIN_NUM = 1;
    VNU_DLY = LDPC_COLUMU_MAX_NUM + LDPC_COLUMU_MIN_NUM + 4;
    % LDPCУ���������1����Ŀ�����ֵ����Сֵ
%     LDPC_ROW_MAX_NUM = max(check_ctrl);%19
%     LDPC_ROW_MIN_NUM = min(check_ctrl);%3
    LDPC_ROW_MAX_NUM = 19;
    LDPC_ROW_MIN_NUM = 3;
    CNU_DLY = LDPC_ROW_MAX_NUM - LDPC_ROW_MIN_NUM + 10;
    SUM_TIME = length(barrel_shifter) + VNU_DLY + 8;

    %% ����У��ڵ���¿�����Ix_c
    %��У��ڵ���¿�����Ix_c�仯Ϊ��ʼ�����������ֹ�������壬��ֹ����-��ʼ����=У��ڵ���¿�����
    %��ʼ����Ixc:__|~~|____________|~~|______________|~~|_______________|~~|_____________
    %��ֹ����Ixd:_______________|~~|______________|~~|_______________|~~|________|~~|____
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

    %% ��������д�뵽coe�ļ���
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
    %% ����ͬZc�Ŀ�����ƴ��

    %% ��������ڵ���¿�����Iy_c
    %�������ڵ���¿�����Iy_c�仯Ϊ��ʼ�����������ֹ�������壬��ֹ����-��ʼ����=�����ڵ���¿�����
    %��ʼ����Iyc:__|~~|____________|~~|______________|~~|_______________|~~|_____________
    %��ֹ����Iyd:_______________|~~|______________|~~|_______________|~~|________|~~|____
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

    %% ��������ڵ㻺���ַ��֯������
    It_a = interweave_addr;     %��ַ��֯������
    %% ��������ڵ㻺������Ͱ����λ��������
    It_c = barrel_shifter;      %��Ͱ����λ������
    It_c_inv = It_c;            %��Ͱ����λ������
    %% ��������ڵ���¶�ȡLLR������
    Iz_c = Iyc;%��ȡLLR����RAMʹ��
    Iz_a1 = zeros(1,length(Iyc));%��ȡLLR����RAM��ַ
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
    %% ƴ�ӱ����ڵ���¿�����Iy_c����ַ����֯������It_a����ַ����֯������It_a_dly����Ͱ����λ������It_c1/It_c3/It_c3����Ͱ����λ������It_c1_dly/It_c3_dly/It_c3_dly
    node_rden   = [  ones(1,len_Iyc),  zeros(1,SUM_TIME - len_Iyc)  ];              %�ڵ㻺���ʹ��
    node_rdaddr = [  It_a           ,  zeros(1,SUM_TIME - len_Iyc)  ];              %�ڵ㻺�����ַ

    barrel_c1   = [zeros(1,2),It_c,zeros(1,SUM_TIME - len_Iyc - 2)];                %��Ͱ����λ�����֣�һ����

    llr_rden    = [zeros(1,4),Iz_c ,zeros(1,SUM_TIME - len_Iyc - 4)];               %LLR�����ʹ��
    llr_rdaddr  = [zeros(1,4),Iz_a1,zeros(1,SUM_TIME - len_Iyc - 4)];               %LLR�������ַ

    llr_dvld    = [zeros(1,5),Iz_c,           zeros(1,SUM_TIME - len_Iyc - 5)];     %I_llr_dvld
    node_dvld   = [zeros(1,5),ones(1,len_Iyc),zeros(1,SUM_TIME - len_Iyc - 5)];     %I_node_dvld
    node_dneg   = [zeros(1,5),Iyd,            zeros(1,SUM_TIME - len_Iyc - 5)];     %I_node_dneg

    barrel_c1_inv = [zeros(1,5+VNU_DLY+0),It_c_inv,       zeros(1,SUM_TIME - len_Iyc - (5+VNU_DLY+0))]; %��Ͱ����λ�����֣�һ����
    node_wren     = [zeros(1,5+VNU_DLY+3),ones(1,len_Iyc),zeros(1,SUM_TIME - len_Iyc - (5+VNU_DLY+3))]; %�ڵ㻺��дʹ��
    node_wraddr   = [zeros(1,5+VNU_DLY+3),It_a,           zeros(1,SUM_TIME - len_Iyc - (5+VNU_DLY+3))]; %�ڵ㻺��д��ַ


    %% ƴ�ӱ����ڵ���¿�����Iy_c����ַ����֯������It_a����ַ����֯������It_a_dly����Ͱ����λ������It_c1/It_c3/It_c3����Ͱ����λ������It_c1_dly/It_c3_dly/It_c3_dly
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
    %% ����ͬZc�Ŀ�����ƴ��

end
    %% �������ִ洢Ϊcoe��ʽ
    file_Ix_c = strcat(Ix_c_name_str,'.coe');
    ptr_Ix_c  = fopen(file_Ix_c,'w');
        fprintf(ptr_Ix_c,'memory_initialization_radix=%d;\n',10);
        fprintf(ptr_Ix_c,'memory_initialization_vector= \n');
        fprintf(ptr_Ix_c,'%d,\n',Ixc_Ixd(1:end-1));
        fprintf(ptr_Ix_c,'%d;\n',Ixc_Ixd(end));
    fclose(ptr_Ix_c);
    %% �������ִ洢Ϊcoe��ʽ
    file_Iyc_Ita_Itc = strcat(Iyc_Ita_Itc_name_str,'.coe');
    buf_Iyc_Ita_Itc = buf_Iyc_Ita_Itc';
    ptr_Iyc_Ita_Itc  = fopen(file_Iyc_Ita_Itc,'w');
        fprintf(ptr_Iyc_Ita_Itc,'memory_initialization_radix=%d;\n',10);
        fprintf(ptr_Iyc_Ita_Itc,'memory_initialization_vector= \n');
        fprintf(ptr_Iyc_Ita_Itc,'%d,\n',buf_Iyc_Ita_Itc(1:end-1));
        fprintf(ptr_Iyc_Ita_Itc,'%d;\n',buf_Iyc_Ita_Itc(end));
    fclose(ptr_Iyc_Ita_Itc);
    display('COE generate finish!')
