%最小和算法，360度并行，1/2码率
function y = ldpc_05_Decoder(llr_data,intera,LDPC_N,LDPC_K, LDPC_P, LDPC_Q,inter_wv_addr,shuffer_ctrl,variable_ctrl,check_ctrl,Zc)

% load('ldpc_1_2_param.mat'); 
% load('ldpc_1_2_param_wbq.mat'); 
% 
% LDPC_N = 16200;  %码长
% LDPC_K = 7200;   %信息长度
% LDPC_P = LDPC_N-LDPC_K;    %校验码长度
% LDPC_Q = (LDPC_N-LDPC_K)/360;
sigma = 0.85;
% sigma = 1;

%%
%开始译码模块

%LLR缓存调整
llr_ram_signal = (reshape(llr_data(1:LDPC_K),Zc,LDPC_K/Zc))';
% llr_ram_check = reshape(llr_data(LDPC_K+1:LDPC_N),LDPC_Q,LDPC_P/LDPC_Q);
llr_ram_check = (reshape(llr_data(LDPC_K+1:LDPC_N),Zc,LDPC_P/Zc))';
llr_ram = [llr_ram_signal;llr_ram_check];

%%
%节点缓存初始化
node_buf = zeros(sum(check_ctrl),Zc);

%%
%迭代
for intera_cnt = intera : -1 : 1 %迭代次数
%%
%变量节点更新/初始化节点缓存

variable_base = 0;

for llr_num = 1:LDPC_N/Zc
    %取数求总和
    colum_sum = llr_ram(llr_num,:); %列求和初始化为llr
    
    for variable_num = 1:variable_ctrl(llr_num)  %variable_ctrl = 8 8 8 8 8 3 3 3......
       buf_data = node_buf(inter_wv_addr(variable_base + variable_num),:); %取出交织地址后的135组数据
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
node_buf(node_buf > 31) = 31;%限幅
node_buf(node_buf < -31) = -31;

%%
%校验节点更新

check_base = 0; %校验节点缓存基地址

for check_num = 1:LDPC_Q
    %取数求最小值和符号
    buf_temp2 = zeros(73,Zc);

    buf_temp1 = node_buf(check_base + 1:check_base + check_ctrl(check_num),:);%取出节点缓存中的一组数据 
    
    sign_all = prod(sign(buf_temp1)); %每行符号相乘，360个得出总符号
    for check_once = 1:check_ctrl(check_num)  %循环 5 6 7 6...
        buf_temp3 = buf_temp1;        
        buf_temp1(buf_temp1 == 0) = 1; %将矩阵中的0替换成1，可正确计算数据符号       
        buf_temp3(check_once,:) = []; %将自己先更新为最大正值，然后更新为其他所有数据的最小值
        
        buf_temp2(check_base + check_once,:) = sign_all.*sign(buf_temp1(check_once,:)).* min(abs(buf_temp3),[],1); %除自己，每行最小值乘以符号
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

