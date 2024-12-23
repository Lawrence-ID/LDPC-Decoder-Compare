% 5g ldpc encoding
% input: s, bit sequence of dimension K * 1
% output: encoded_bits, bit sequence
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT

%函数输入：输入待编码比特序列s，输入基图索引
%函数输出：输出编码比特s、校验矩阵H、ZC数值、编码后的二进制序列(s,w)

function [encoded_bits, H, Z_c, encoded_bits_original] = ldpc_encode(s, base_graph_index,kb)

K = length(s);  %求出待编码序列的长度

encoded_bits = zeros(3*K, 1);%生成一个列向量，行数为3K
%BG1的大小为46*68

if base_graph_index == 1    %如果基图选择1
    a = 4;                  %BG1H的核心矩阵部分为4*22
    b = 22;                 %BG1的可编码的比特为22Zc
    c = 26;
    d = 42;
    e = 46;
    z = K/b;                %Zc的数值，标准的计算方法，这个数值需要处理，计算出的Zc数值可能不是标准中所定义的，需要进行扩充
    Z_c = z;
    N = 66 * Z_c;
    z = K/b;
    set_index = lifting_size_table_lookup(z);%主要任务做ILS索引，在标准中定义了八种ILS索引，所以实际的BG1共有八种形式，这个函数的实现实现方法比较神奇，但是缺点是单单支持标准中所定义的几种Z值的ILS索引，当码块的数值不符合标准的数值怎么处理
    load parity_check_matrices_protocol_1    %文件中存储了BG1所有基图格式，索引的正确数值范围为1~8
    BG = parity_check_matrices_protocol_1(:, :, set_index); %#ok<NODEF>   选出BG1的基图
%BG1的大小为42*52

elseif base_graph_index == 2%          
    a = 4;                  %BG2H的核心矩阵部分为4*22
    b = 10;                 %BG2的可编码的比特为10Zc
    c = 14;
    d = 38;
    e = 42;
    z = K/b;                %Zc的数值为
    Z_c = z;
    N = 50 * Z_c;
    set_index = lifting_size_table_lookup(z); %主要任务做ILS索引，在标准中定义了八种ILS索引，所以实际的BG2共有八种形式   协议38 212 19页
    load parity_check_matrices_protocol_2     %文件中存储了BG2所有基图格式，索引的正确数值范围为1~8
    BG = parity_check_matrices_protocol_2(:, :, set_index); %#ok<NODEF>     选出BG2的基图
else
  error('wrong base graph index in ldpc encoding.');
end

BG(BG ~= -1) = mod(BG(BG ~= -1), Z_c); %骚操作，是不是看不懂，其实很简单，就是把选好的基图对Zc取余，计算出BG所对应的循环移位矩阵
%函数功能打掉BG前两列，并且把-1元素转换为0元素---意思是BG的前两列打掉不传输---但是实际编码还会用到


% set_index = lifting_size_table_lookup(Z_c);
% BG = parity_check_matrices_protocol(:, :, set_index);
%下面是把BG矩阵进行分块处理
%具体的划分形式见本文件夹图片
A_prime = BG(1:a, 1:b);             %矩阵A--核心矩阵部分    a*b
B_prime = BG(1:a, (b+1):c);         %矩阵B--双对角结构矩阵      a*a
C_prime = BG((a+1):e, 1:b);         %矩阵C--扩展矩阵（a single parity-check）
D_prime = BG((a+1):e, (b+1):c);     %矩阵D--扩展矩阵（a single parity-check）  C+D

if base_graph_index == 2 %在bg2时修改b值
    if(kb == 10)        
        b = 10;
        s = s;
        K = K;
    elseif(kb == 9)
        b = 9;
        s = s(1:9*Z_c);
        K = K - 1 * Z_c;
    elseif(kb == 8)
        b = 8;
        s = s(1:8*Z_c);
        K = K - 2 * Z_c;
    elseif(kb == 6)
        b = 6;
        s = s(1:6*Z_c);
        K = K - 4 * Z_c;
    else
        error('wrong kb in ldpc encoding.');
    end
end

%根据输入的kb取值，修正矩阵A和C的取值范围
A_prime = A_prime(:, 1:b);             %矩阵A--核心矩阵部分    a*b
C_prime = C_prime(:, 1:b);         %矩阵C--扩展矩阵（a single parity-check）

encoded_bits(1:K-2*Z_c) = s(2*Z_c+1:end);
% for k = (2*Z_c):(K-1)
%   if s(k+1) ~= -1
%     encoded_bits(k-2*Z_c+1) = s(k+1);
%   else
%     s(k+1) = 0;
%     encoded_bits(k-2*Z_c+1) = -1;
%   end    
% end

z = Z_c;
%spalloc函数功能：为稀疏矩阵分配空间--spalloc(10,10,20)---非零元素最多为20
%A的大小为--a*z×b*z-
%size(A_prime)返回a*z和b*z
%ones(size(A_prime))---产生一个a*z×b*z大小的单位矩阵
%为什么这么找出这个矩阵中非零最大值??
A = spalloc(a*z, b*z, nnz(A_prime + ones(size(A_prime))));
for row_index = 1:a     %行索引  1:a
    for column_index = 1:b      %列索引     1:b
        if A_prime(row_index, column_index) ~= -1  %if(矩阵内部的数值不等于-1)  begin---下面的过程是将A矩阵中的非零元素--根据标准定义的---通过循环移位---得到扩展后的矩阵
            A((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(A_prime(row_index, column_index), z)+1):z, 1:mod(A_prime(row_index, column_index), z)], ones(1, z), z, z);
        end
    end
end
%矩阵B本身就是4*4大小
B = spalloc(a*z, a*z, nnz(B_prime + ones(size(B_prime))));

for row_index = 1:a
    for column_index = 1:a
        if B_prime(row_index, column_index) ~= -1
            B((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(B_prime(row_index, column_index), z)+1):z, 1:mod(B_prime(row_index, column_index), z)], ones(1, z), z, z);
        end
    end
end
% d=e-a  剩下的行数
C = spalloc(d*z, b*z, nnz(C_prime + ones(size(C_prime))));

for row_index = 1:d
    for column_index = 1:b
        if C_prime(row_index, column_index) ~= -1
            C((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(C_prime(row_index, column_index), z)+1):z, 1:mod(C_prime(row_index, column_index), z)], ones(1, z), z, z);
        end
    end
end

D = spalloc(d*z, a*z, nnz(D_prime + ones(size(D_prime))));

for row_index = 1:d
    for column_index = 1:a
        if D_prime(row_index, column_index) ~= -1
            D((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(D_prime(row_index, column_index), z)+1):z, 1:mod(D_prime(row_index, column_index), z)], ones(1, z), z, z);
        else
            D((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = spalloc(z, z, 0);
        end
    end
end

B_inv = spalloc(a*z, a*z, 20*z);

%speye函数功能：返回一个az*az的对角矩阵，求的是矩阵的逆变换B的逆变换
%生成一块矩阵  见下面文件
%每个基图对应的4*4具有双对角结构的矩阵的元素均不相同，这部分需要特殊处理
if (base_graph_index == 1) && (set_index ~= 7)    %如果选用BG1---set_index！=7---BG1基图的双对角结构部分有两种，这里分开计算
    B_inv(1:z, 1:z)             = speye(z);       %B_inv矩阵的z*z位置为单位矩阵
    B_inv(1:z, 1+z:2*z)         = speye(z);       %两个z*z的对角矩阵组合  大小为z*2z
    B_inv(1:z, 1+2*z:3*z)       = speye(z);       %三个z*z的对角矩阵组合  大小为z*3z
    B_inv(1:z, 1+3*z:4*z)       = speye(z);       %四个z*z的对角矩阵组合  大小为z*4z
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
%单位阵的基础上向右移动z-1位
elseif (base_graph_index == 2) && ((set_index ~= 4) && (set_index ~= 8))
    B_inv(1:z, 1:z)             = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+z:2*z)         = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+2*z:3*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+3*z:4*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
elseif (base_graph_index == 1) && (z == 208)
    B_inv(1:z, 1:z)             = sparse(circshift(eye(208), 105));
    B_inv(1:z, 1+z:2*z)         = sparse(circshift(eye(208), 105));
    B_inv(1:z, 1+2*z:3*z)       = sparse(circshift(eye(208), 105));
    B_inv(1:z, 1+3*z:4*z)       = sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1:z)         = speye(208) + sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1:z)       = sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1+z:2*z)   = sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1+2*z:3*z) = speye(208) + sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1+3*z:4*z) = speye(208) + sparse(circshift(eye(208), 105));
    B_inv(1+3*z:4*z, 1:z)       = sparse(circshift(eye(208), 105)); 
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(circshift(eye(208), 105)); 
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(circshift(eye(208), 105)); 
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(208) + sparse(circshift(eye(208), 105));
elseif (base_graph_index == 1) && ((z ~= 208) && (set_index == 7))
    B_inv(1:z, 1:z)             = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+z:2*z)         = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+2*z:3*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+3*z:4*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
elseif (base_graph_index == 2) && ((set_index == 4) || (set_index == 8))    
    B_inv(1:z, 1:z)             = speye(z);
    B_inv(1:z, 1+z:2*z)         = speye(z);
    B_inv(1:z, 1+2*z:3*z)       = speye(z);
    B_inv(1:z, 1+3*z:4*z)       = speye(z);
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);    
end

s = s(:);
%公式---近似下三角矩阵编码的方法   B_inv为B的逆矩阵---
p_1 = mod(B_inv * (A * s), 2);
p_2 = mod(C * s + D * p_1, 2);
w = [p_1; p_2];

%把编码的比特补全校验位
% for k = K:(N+2*Z_c-1)
%    encoded_bits(k-2*Z_c+1) = w(k-K+1);   
% end    
encoded_bits(K-2*Z_c+1:K-2*Z_c+length(w)) = w;
%生成H矩阵  前提已经打去了前两列的BG---BG1 46*66  BG2---42*52，生成的H矩阵
%spalloc(a*z, d*z, 0)生成0阵        speye(d*z)生成对角矩阵
H = [A, B, spalloc(a*z, d*z, 0); C, D, speye(d*z)];

%原始的编码比特，BG对应的前两列保留的情况
encoded_bits_original = [s; w]; 

% mod(H * [s; p_1; p_2]) = 0

clear A
clear B
clear C
clear D
clear B_inv
clear A_prime
clear B_prime
clear C_prime
clear D_prime

end

%前两个矩阵的测试  Z=4，生成结果如下图所示：
%1     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0
%0     1     0     0     0     1     0     0     0     1     0     0     0     1     0     0
%0     0     1     0     0     0     1     0     0     0     1     0     0     0     1     0
%0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     1
%1     1     0     0     0     1     0     0     0     1     0     0     0     1     0     0
%0     1     1     0     0     0     1     0     0     0     1     0     0     0     1     0
%0     0     1     1     0     0     0     1     0     0     0     1     0     0     0     1
%1     0     0     1     1     0     0     0     1     0     0     0     1     0     0     0
%0     1     0     0     0     1     0     0     1     1     0     0     1     1     0     0
%0     0     1     0     0     0     1     0     0     1     1     0     0     1     1     0
%0     0     0     1     0     0     0     1     0     0     1     1     0     0     1     1
%1     0     0     0     1     0     0     0     1     0     0     1     1     0     0     1
%0     1     0     0     0     1     0     0     0     1     0     0     1     1     0     0
%0     0     1     0     0     0     1     0     0     0     1     0     0     1     1     0
%0     0     0     1     0     0     0     1     0     0     0     1     0     0     1     1
%1     0     0     0     1     0     0     0     1     0     0     0     1     0     0     1


%向右移动3位
%0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     1
%1     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0
%0     1     0     0     0     1     0     0     0     1     0     0     0     1     0     0
%0     0     1     0     0     0     1     0     0     0     1     0     0     0     1     0
%1     0     0     1     0     0     0     1     0     0     0     1     0     0     0     1
%1     1     0     0     1     0     0     0     1     0     0     0     1     0     0     0
%0     1     1     0     0     1     0     0     0     1     0     0     0     1     0     0
%0     0     1     1     0     0     1     0     0     0     1     0     0     0     1     0
%1     0     0     1     1     0     0     1     0     0     0     1     0     0     0     1
%1     1     0     0     1     1     0     0     1     0     0     0     1     0     0     0
%0     1     1     0     0     1     1     0     0     1     0     0     0     1     0     0
%0     0     1     1     0     0     1     1     0     0     1     0     0     0     1     0
%0     0     0     1     0     0     0     1     0     0     0     1     1     0     0     1
%1     0     0     0     1     0     0     0     1     0     0     0     1     1     0     0
%0     1     0     0     0     1     0     0     0     1     0     0     0     1     1     0
%0     0     1     0     0     0     1     0     0     0     1     0     0     0     1     1
