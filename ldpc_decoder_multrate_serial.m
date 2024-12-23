function [y, llr_final] = ldpc_decoder_multrate_serial(llr_data,llr_width,intera,LDPC_N,LDPC_K, LDPC_P, LDPC_Q,inter_wv_addr,shuffer_ctrl,variable_ctrl,check_ctrl,Zc)

    % 描述校验矩阵的分割方法，以BG2、Zc=104为例
    % 1. 根据算法特征，将校验矩阵按照行进行分组，每间隔104行，拿出来一行放一组，所以一共可以分成104组，每组共有42行
    % 2. 每一组中的1的个数是相同的，因为每个Zc块是单位矩阵循环移位得出的。这里1的个数是197，所以可以使用104×197矩阵存储LLR
    
    % 描述非零元素的相关索引方法，以BG2、Zc=104为例
    % 1. H_T矩阵是对H_matrix进行行变换的过程，H_matrix的行每隔104个抽出一行分成一组，可以分成104组，每组共有42行
    % 2. 找出H_T矩阵每行的非零元素索引，存储在buf_a中
    % 3. 将buf_a中的每42行弄成一行，存储在buf_b中
    % 4. 找出buf_b中的非零元素的索引，存储在buf_c中
    
    % 描述变量节点的更新方法，以BG2、Zc=104为例
    % 1. 变量节点的更新是以列为单位进行更新的，首先需要计算出第一组中的每一列的非零元素的个数：variable_ctrl
    % 2. 找到非零元素的个数后，还需要记录按照列排列后的每一个非零元素的索引：inter_wv_addr（1-197）
    % 3. 还需要找到移位系数：shuffer_ctrl
    
    % 产生参数
        % variable_ctrl：把具体的索引数值以Zc为区间进行分组，variable_ctrl为每组索引值的数量。
        % interweave_addr：把buf_d数值先按照Zc进行分区，将分区后的索引以此存在。
        % barrel_shifter：将index对Zc取模，如果计算结果为-1，shuffer_ctrl_s(base+jj) = Zc+temp_shuffer;
        
    sigma = 1;
    
    %%
    %开始译码模块
    
    %LLR缓存调整
    llr_ram_signal = (reshape(llr_data(1:LDPC_K),Zc,LDPC_K/Zc))';
    % llr_ram_check = reshape(llr_data(LDPC_K+1:LDPC_N),LDPC_Q,LDPC_P/LDPC_Q);
    llr_ram_check = (reshape(llr_data(LDPC_K+1:LDPC_N),Zc,LDPC_P/Zc))';
    llr_ram = [llr_ram_signal;llr_ram_check];
    
    %%
    %节点缓存初始化，存储检验矩阵非零元素
    node_buf = zeros(sum(check_ctrl),Zc);
    
    %%
    %迭代
    for intera_cnt = intera : -1 : 1 %迭代次数
    %%
    %变量节点更新/初始化节点缓存
    
    variable_base = 0;

    maxMagnitude = 2^(llr_width-1);
    
    for llr_num = 1:LDPC_N/Zc
        % 取数求总和
        % 注意这部分，每取出104，对应校验矩阵的前104列。
        colum_sum = llr_ram(llr_num,:); %列求和初始化为llr，将LLR的第一行取出来，准备变量节点更新
        
        for variable_num = 1:variable_ctrl(llr_num)  % 校验矩阵中每一列非零元素的个数
           buf_data = node_buf(inter_wv_addr(variable_base + variable_num),:); %将node_buf的第variable_num行取出来
           buf_data_shift = circshift(buf_data,[0 shuffer_ctrl(variable_base + variable_num)]);%将数据正循环移位(后)
           colum_sum = buf_data_shift + colum_sum;%
        end
        
        %存储求和数据
        for variable_num = 1:variable_ctrl(llr_num)  %variable_ctrl = 8 8 8 8 8 3 3 3......      
           %取原缓存数据
           buf_data  = node_buf(inter_wv_addr(variable_base + variable_num),:); %取出交织地址后的135组数据
           buf_data_shift = circshift(buf_data,[0 shuffer_ctrl(variable_base + variable_num)]);%将数据正循环移位(后)
           %用总和减原缓存数据
           buf_data_save = colum_sum - buf_data_shift;  
           %存储到原缓存位置
           node_buf(inter_wv_addr(variable_base + variable_num),:) = circshift(buf_data_save,[0 -shuffer_ctrl(variable_base + variable_num)]);
        end
        
        variable_base = variable_base + variable_ctrl(llr_num);%基地址变化
    end
    
    node_buf_vnu1 = node_buf';
    %%
    %节点缓存 定点
    % node_buf = fi(node_buf,1,10,0);
    % node_buf = node_buf.data;

    % 去零、限幅
    node_buf(node_buf == 0) = 1;%去零
    node_buf(node_buf > maxMagnitude - 1) = maxMagnitude - 1;%限幅
    node_buf(node_buf < -maxMagnitude) = -maxMagnitude;
    
    %%
    %校验节点更新
    
    check_base = 0; %校验节点缓存基地址
    
    for check_num = 1:LDPC_Q
        %取数求最小值和符号
        buf_temp2 = zeros(316,Zc);
    
        buf_temp1 = node_buf(check_base + 1:check_base + check_ctrl(check_num),:);%取出节点缓存中的一组数据 
        
        sign_all = prod(signStrict(buf_temp1)); %每行符号相乘，360个得出总符号
        for check_once = 1:check_ctrl(check_num)  %循环 5 6 7 6...
            buf_temp3 = buf_temp1;        
            buf_temp1(buf_temp1 == 0) = 1; %将矩阵中的0替换成1，可正确计算数据符号       
            buf_temp3(check_once,:) = []; %将自己先更新为最大正值，然后更新为其他所有数据的最小值
            
            buf_temp2(check_base + check_once,:) = sign_all.*signStrict(buf_temp1(check_once,:)).* min(abs(buf_temp3),[],1); %除自己，每行最小值乘以符号
        end
        
        %整组数据进行更新
        node_buf([check_base + 1:check_base + check_ctrl(check_num)],:) = buf_temp2([check_base + 1:check_base + check_ctrl(check_num)],:); 
     
        check_base = check_base + check_ctrl(check_num);%基地址变化
    end
    
    %控制系数
    % if(intera_cnt ~= intera)
        node_buf = sigma.*node_buf;
    % end
    
    node_buf_cnu1 = node_buf';
    end
    %%
    %译码输出
    variable_base = 0;
    decode_data = zeros(LDPC_K/Zc,Zc);
    for llr_num = 1:LDPC_K/Zc
        %取数求总和
        colum_sum = llr_ram(llr_num,:); %列求和初始化为llr
        
        for variable_num = 1:variable_ctrl(llr_num)  %variable_ctrl = 8 8 8 8 8 3 3 3......
           buf_data = node_buf(inter_wv_addr(variable_base + variable_num),:); %取出交织地址后的135组数据
           buf_data_shift = circshift(buf_data,[0 shuffer_ctrl(variable_base + variable_num)]);%将数据正循环移位(后)
           colum_sum = buf_data_shift + colum_sum;
        end
        
        decode_data(llr_num,:) = colum_sum; %求和译码，每次译360个数，存储20组
        
        variable_base = variable_base + variable_ctrl(llr_num);%基地址变化
    end
    
       serial_code = reshape(decode_data',1,LDPC_K);
       llr_final = serial_code;
       
       %数据译码判决
       decode_judge = zeros(1,LDPC_K);
       for i = 1:LDPC_K
            if(serial_code(i) >= 0)
                decode_judge(i) = 0;  
            else
                decode_judge(i) = 1; 
            end
       end
    
       y = decode_judge;
    
    