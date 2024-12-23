%��С���㷨��360�Ȳ��У�1/2����
function y = ldpc_05_Decoder(llr_data,intera,LDPC_N,LDPC_K, LDPC_P, LDPC_Q,inter_wv_addr,shuffer_ctrl,variable_ctrl,check_ctrl,Zc)

% load('ldpc_1_2_param.mat'); 
% load('ldpc_1_2_param_wbq.mat'); 
% 
% LDPC_N = 16200;  %�볤
% LDPC_K = 7200;   %��Ϣ����
% LDPC_P = LDPC_N-LDPC_K;    %У���볤��
% LDPC_Q = (LDPC_N-LDPC_K)/360;
sigma = 0.85;
% sigma = 1;

%%
%��ʼ����ģ��

%LLR�������
llr_ram_signal = (reshape(llr_data(1:LDPC_K),Zc,LDPC_K/Zc))';
% llr_ram_check = reshape(llr_data(LDPC_K+1:LDPC_N),LDPC_Q,LDPC_P/LDPC_Q);
llr_ram_check = (reshape(llr_data(LDPC_K+1:LDPC_N),Zc,LDPC_P/Zc))';
llr_ram = [llr_ram_signal;llr_ram_check];

%%
%�ڵ㻺���ʼ��
node_buf = zeros(sum(check_ctrl),Zc);

%%
%����
for intera_cnt = intera : -1 : 1 %��������
%%
%�����ڵ����/��ʼ���ڵ㻺��

variable_base = 0;

for llr_num = 1:LDPC_N/Zc
    %ȡ�����ܺ�
    colum_sum = llr_ram(llr_num,:); %����ͳ�ʼ��Ϊllr
    
    for variable_num = 1:variable_ctrl(llr_num)  %variable_ctrl = 8 8 8 8 8 3 3 3......
       buf_data = node_buf(inter_wv_addr(variable_base + variable_num),:); %ȡ����֯��ַ���135������
       buf_data_shift = circshift(buf_data,[0 shuffer_ctrl(variable_base + variable_num)]);%��������ѭ����λ(��)
       colum_sum = buf_data_shift + colum_sum;%
    end
    
    %�洢�������
    for variable_num = 1:variable_ctrl(llr_num)  %variable_ctrl = 8 8 8 8 8 3 3 3......      
       %ȡԭ��������
       buf_data  = node_buf(inter_wv_addr(variable_base + variable_num),:); %ȡ����֯��ַ���135������
       buf_data_shift = circshift(buf_data,[0 shuffer_ctrl(variable_base + variable_num)]);%��������ѭ����λ(��)
       %���ܺͼ�ԭ��������
       buf_data_save = colum_sum - buf_data_shift;  
       %�洢��ԭ����λ��
       node_buf(inter_wv_addr(variable_base + variable_num),:) = circshift(buf_data_save,[0 -shuffer_ctrl(variable_base + variable_num)]);
    end
    
    variable_base = variable_base + variable_ctrl(llr_num);%����ַ�仯
end

node_buf_vnu1 = node_buf';
%%
%�ڵ㻺�� ����
% node_buf = fi(node_buf,1,10,0);
% node_buf = node_buf.data;

% ȥ�㡢�޷�
node_buf(node_buf == 0) = 1;%ȥ��
node_buf(node_buf > 31) = 31;%�޷�
node_buf(node_buf < -31) = -31;

%%
%У��ڵ����

check_base = 0; %У��ڵ㻺�����ַ

for check_num = 1:LDPC_Q
    %ȡ������Сֵ�ͷ���
    buf_temp2 = zeros(73,Zc);

    buf_temp1 = node_buf(check_base + 1:check_base + check_ctrl(check_num),:);%ȡ���ڵ㻺���е�һ������ 
    
    sign_all = prod(sign(buf_temp1)); %ÿ�з�����ˣ�360���ó��ܷ���
    for check_once = 1:check_ctrl(check_num)  %ѭ�� 5 6 7 6...
        buf_temp3 = buf_temp1;        
        buf_temp1(buf_temp1 == 0) = 1; %�������е�0�滻��1������ȷ�������ݷ���       
        buf_temp3(check_once,:) = []; %���Լ��ȸ���Ϊ�����ֵ��Ȼ�����Ϊ�����������ݵ���Сֵ
        
        buf_temp2(check_base + check_once,:) = sign_all.*sign(buf_temp1(check_once,:)).* min(abs(buf_temp3),[],1); %���Լ���ÿ����Сֵ���Է���
    end
    
    %�������ݽ��и���
    node_buf([check_base + 1:check_base + check_ctrl(check_num)],:) = buf_temp2([check_base + 1:check_base + check_ctrl(check_num)],:); 
 
    check_base = check_base + check_ctrl(check_num);%����ַ�仯
end

%����ϵ��
% if(intera_cnt ~= intera)
    node_buf = sigma.*node_buf;
% end

node_buf_cnu1 = node_buf';
end
%%
%�������
variable_base = 0;
decode_data = zeros(LDPC_K/Zc,Zc);
for llr_num = 1:LDPC_K/Zc
    %ȡ�����ܺ�
    colum_sum = llr_ram(llr_num,:); %����ͳ�ʼ��Ϊllr
    
    for variable_num = 1:variable_ctrl(llr_num)  %variable_ctrl = 8 8 8 8 8 3 3 3......
       buf_data = node_buf(inter_wv_addr(variable_base + variable_num),:); %ȡ����֯��ַ���135������
       buf_data_shift = circshift(buf_data,[0 shuffer_ctrl(variable_base + variable_num)]);%��������ѭ����λ(��)
       colum_sum = buf_data_shift + colum_sum;
    end
    
    decode_data(llr_num,:) = colum_sum; %������룬ÿ����360�������洢20��
    
    variable_base = variable_base + variable_ctrl(llr_num);%����ַ�仯
end

   serial_code = reshape(decode_data',1,LDPC_K);
   
   %���������о�
   decode_judge = zeros(1,LDPC_K);
   for i = 1:LDPC_K
        if(serial_code(i) >= 0)
            decode_judge(i) = 0;  
        else
            decode_judge(i) = 1; 
        end
   end

   y = decode_judge;

